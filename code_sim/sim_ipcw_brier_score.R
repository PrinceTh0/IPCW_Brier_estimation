
brier_score <- function(matrix_delta, matrix_obs_times, matrix_ref_times, matrix_G_hat_event_t_min_1, G_hat_t, S_hat, n_test, max_weight = 100, brier_general = T){
  
  I_failed <- ifelse((matrix_delta == 1) & (matrix_obs_times <= matrix_ref_times), 1, 0)
  I_survived <- ifelse((matrix_obs_times > matrix_ref_times), 1, 0)
  
  w_failed <- I_failed/matrix_G_hat_event_t_min_1
  w_survived <- I_survived/G_hat_t
  
  #w_failed <- ifelse(w_failed > max_weight, max_weight, w_failed)
  #w_survived <- ifelse(w_survived > max_weight, max_weight, w_survived)
  
  weights <- w_failed + w_survived
  
  n_tilde <- colSums(weights)
  
  failed <- ((S_hat)**2) * w_failed
  survived <- ((1-S_hat)**2) * w_survived
  
  if(brier_general == T){
    brier <- (1/n_tilde)*colSums(failed + survived)
  } else{
    brier <- (1/n_test)*colSums(failed + survived)
  }
  
  return(brier)
}
