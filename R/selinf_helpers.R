# Approximate MLE for data splitting, written by Jonathan Taylor
approximate_mle_inference = function(training_proportion,
                                     training_betahat,
                                     selected_beta_refit,
                                     selected_signs,
                                     selected_hessian,
                                     selected_feature_weights,
                                     level=0.9) {

    nselect = nrow(selected_hessian)
    pi_s = training_proportion
    ratio = (1 - pi_s) / pi_s

    target_cov = chol2inv(chol(selected_hessian))
    #cond_precision = selected_hessian / ratio
    #cond_cov = target_cov * ratio
    Signs = diag(selected_signs)                                          
    #cond_cov = Signs %*% cond_cov %*% Signs
    prec_opt = (Signs %*% selected_hessian %*% Signs)/ratio

    logdens_linear = Signs %*% target_cov 
    cond_mean = (selected_beta_refit * selected_signs - logdens_linear %*% (
                 selected_feature_weights *
                 selected_signs))
    linear_part = -diag(rep(1, nselect))
    offset = rep(0, nselect)

    target_score_cov = -diag(rep(1, nselect))
    observed_target = selected_beta_refit
    
    result = selective_MLE(observed_target, 
                           target_cov,
                           target_score_cov, 
                           training_betahat * selected_signs,
                           cond_mean,
                           prec_opt,
                           logdens_linear,
                           linear_part,
                           offset,
                           selected_hessian,
                           level=level)

    return(result)
}


selective_MLE = function(observed_target, 
                         target_cov,
                         target_score_cov, 
                         feasible_point,
                         cond_mean,
                         prec_opt,
                         logdens_linear,
                         linear_part,
                         offset,
                         prec_target,
                         level=0.9,
                         step=1,
                         max_iter=1000,
                         min_iter=200,
                         tol=1.e-12) {

    #prec_target = chol2inv(chol(target_cov))
    target_lin = - logdens_linear %*% t(target_score_cov) %*% prec_target
    target_offset = cond_mean - target_lin %*% observed_target

    #prec_opt = chol2inv(chol(cond_cov))
    conjugate_arg = prec_opt %*% cond_mean

    solve_result = solve_barrier_affine(conjugate_arg,
                                        prec_opt,
                                        feasible_point,
                                        linear_part,
                                        offset,
                                        step=step,
                                        max_iter=max_iter,
                                        min_iter=min_iter,
                                        tol=tol) 

    soln = solve_result$opt_variable
    hess = solve_result$hess

    final_estimator = observed_target + target_cov %*% t(target_lin) %*% prec_opt  %*% (cond_mean - soln)
    ind_unbiased_estimator = observed_target + target_cov %*% t(target_lin) %*% prec_opt %*%  (cond_mean - feasible_point)

    L = t(target_lin) %*% prec_opt
    observed_info_natural = prec_target + L %*% target_lin - L %*% hess %*% t(L)
    observed_info_mean = target_cov %*% observed_info_natural %*% target_cov

    Z_scores = final_estimator / sqrt(diag(observed_info_mean))
    pvalues = pnorm(Z_scores)
    pvalues = 2 * pmin(pvalues,  1 - pvalues)

    alpha = 1 - level
    q = qnorm(1 - alpha / 2.)
    se_ = sqrt(diag(observed_info_mean))

    result = data.frame(MLE=final_estimator,
                        SE=se_,
                        Zvalue=Z_scores,
                        pvalue=pvalues,
                        lower_confidence=final_estimator - q * se_,
                        upper_confidence=final_estimator + q * se_,
                        unbiased=ind_unbiased_estimator)
    return(list(summary=result, vcov=observed_info_mean))
}                             

solve_barrier_affine = function(conjugate_arg,
                                precision,
                                feasible_point,
                                linear_term,
                                offset,
                                step=1,
                                max_iter=1000,
                                min_iter=200,
                                tol=1.e-12) {

    gradient = rep(0, length(conjugate_arg))
    opt_variable = 1 * feasible_point
    opt_proposed = opt_variable  * 1 # force a copy
    affine_term = 0 * offset
    A = linear_term

    scaling = sqrt(diag(A %*% precision %*% t(A))) # linear term is -I for LASSO

    result = selectiveInference:::solve_barrier_(conjugate_arg,
                                                 precision,	
                                                 feasible_point,
                                                 max_iter,
                                                 min_iter,
                                                 tol,
                                                 step)

    opt_variable = result$soln
    final_affine = offset + opt_variable # linear_term is -I for LASSO
    diag_barrier = - 1 / (final_affine + scaling)^2 + 1 / final_affine^2
    if(length(opt_variable) > 1) {
       hess = chol2inv(chol(precision + diag(diag_barrier)))
    } else {
       hess = 1 / (precision + diag_barrier)
       hess = matrix(hess, 1, 1)
    }
    return(list(opt_variable=opt_variable,
                hess=hess))
}



rev_cumsum = function(v){
  rev(cumsum(rev(v)))
}

# fast for large n, actually always faster.
# also assume t1 <= t2 <= t3 ...
fast_cox_fisher = function(f, x, time, c){
  f = scale(f,TRUE,FALSE)
  r = rank(time, ties.method="min")
  rm =rank(time, ties.method="max")
  theta = exp(f)[,1]
  rskden=rev_cumsum(theta)[r]
  a = (cumsum(c/rskden)[rm])*theta
  first = sweep(x, 1, sqrt(a), FUN="*")
  first = crossprod(first)
  cum_matrix = sweep(x, 1, theta, FUN="*")
  cum_matrix = apply(cum_matrix, 2, rev_cumsum)[r,]
  second = sweep(cum_matrix, 1, c/rskden, FUN="*")
  second=crossprod(second)
  first - second
}

# A wrapper of the previous function
cox_fisher = function(f, x, time, c){
  o = order(time)
  f = f[o]
  x = x[o,]
  c = c[o]
  time = time[o]
  fast_cox_fisher(f, x, time, c)
  
}
