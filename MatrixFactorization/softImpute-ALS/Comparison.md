Here are the comparison with the software of Zhongyuan Lyu in GitHub and my suggestions.

## Comparison

### Visualization

Lyu is more proficient in the using of package `ggplot2` than I am. What deserves my attention most is the convenience of factor variables in plotting. For example:
```r
res_df <- data.frame(x=c(als_time[m:length(als_time)],
                         sI_time[m:length(sI_time)],
                         sI_als_time[m:length(sI_als_time)]),
                     y=c(als_res$rel_obj[m:length(als_time)],
                         sI_res$rel_obj[m:length(sI_time)],
                         sI_als_res$rel_obj[m:length(sI_als_time)]),
                     method=factor(c(rep(method_name[1],length(als_time)-m+1),
                                     rep(method_name[2],length(sI_time)-m+1),
                                     rep(method_name[3],length(sI_als_time)-m+1))))

# When using function ggplot():
aes(x=x,y=y,color=method)
geom_point(aes(shape = method))
```
It can easily distinguish the results of different methods and show them on one picture.

The output in Lyu's main R code is based on the extensive use of the function `cat()`. For example:
```r
cat("ALS converges in", iter, "iterations")
cat("\n")
```

All in all, I think Lyu's output is more intuitive and detailed than mine. This is the most valuable and important thing that I learned from Lyu's software.

### Convergence condition

In Lyu's software, he calculated the value of the objective function after each update and compared them in log space as the convergence condition:
```r
ALS_obj <- function(X,A,B,lambda){
  0.5*sum((X-A%*%t(B))^2,na.rm = TRUE)+0.5*lambda*(norm(A,"F")^2+norm(B,"F")^2)
}
obj_prev <- ALS_obj(X_obs,A_int,B_int,lambda)
obj <- ALS_obj(X_obs,A,B,lambda)
rel_obj[iter] <- log((obj_prev)/obj)
```
Please notice that Lyu divided the initial value by each term after update.

I compared the F-norm of the solution as the convergence condition.

### Time

The time spent on each update is important data, but I didn't record it.
```r
start_time <- proc.time()[[3]]
time_vec[iter] <- proc.time()[[3]]-start_time
```

### Matrix multiplication

In the update process, we need to do matrix multiplication: `U * D * V`, where `D` is a diagonal matrix. Lyu used `sweep()` in his software:
```r
S <- soft_thred(svd_res$d,lambda)
Mhat <- sweep(svd_res$u, MARGIN = 2, S, "*") %*% t(svd_res$v)
```
I directly used the matrix multiplication function `%*%` in R like this:
```r
Mhat <- U %*% diag(S) %*% t(V)
```

## Suggestion

### `function: ALS` lines 27 to 38

I think the using of for loop will slow down Lyu's software:
```r
for (i in 1:m){
    nz_cols <- c(1:n)[ind_mat[i,]]
    B_submat <- B[nz_cols,]
    BtB <- t(B_submat)%*%B_submat
    A[i,] <- solve(BtB+lambda*diag(r))%*%(t(B_submat)%*%X_obs[i,nz_cols])
}
for (j in 1:n){
    nz_rows <- c(1:m)[ind_mat[,j]]
    A_submat <- A[nz_rows,]
    AtA <- t(A_submat)%*%A_submat
    B[j,] <- solve(AtA+lambda*diag(r))%*%(t(A_submat)%*%X_obs[nz_rows,j]) 
}
```
If it was me, it would be handled like this:
```r
# m, n, X_obs, ind_mat, lambda and r have been set up
UpdateA <- function(i, B){
    nz_cols <- c(1 : n)[ind_mat[i, ]]
    B_submat <- B[nz_cols, ]
    BtB <- t(B_submat) %*% B_submat
    A_i <- solve(BtB + lambda * diag(r)) %*% (t(B_submat) %*%  X_obs[i, nz_cols])
}
UpdateB <- function(j, A){
    nz_rows <- c(1 : m)[ind_mat[ , j]]
    A_submat <- A[nz_rows, ]
    AtA <- t(A_submat) %*% A_submat
    B_j <- solve(AtA + lambda * diag(r)) %*% (t(A_submat) %*%  X_obs[nz_rows, j])
}
A <- t(sapply(1 : m, FUN = UpdateA, B = B))
B <- t(sapply(1 : n, FUN = UpdateB, A = A))
```
In fact, there are some other for loops that can be improved in Lyu's software:
```r
for (i in 1:length(lambda_list)){
  # als_res <- ALS(X_obs, r=r, itermax=itermax, lambda=lambda_list[i], tol=tol)
  # sI_res <- softImpute(X_obs, itermax=itermax, lambda=lambda_list[i], tol=tol)
  # X_als <- als_res$A_est%*%t(als_res$B_est)
  # X_sI <- sI_res$Mhat
  sI_als_res <- softImpute_ALS(X_obs, r=r, itermax=itermax, lambda=lambda_list[i], tol=tol)
  X_sI_als <- sI_als_res$U%*%diag(sI_als_res$Dsigmalambda)%*%t(sI_als_res$V)
  ## Training error
  mse_mat[i, 1] <- mse(!na_mat, X_sI_als, X_true)
  ## Test error
  mse_mat[i, 2] <- mse(na_mat, X_sI_als, ABt_true)
  cat("iteration", i,"\n")
}
```