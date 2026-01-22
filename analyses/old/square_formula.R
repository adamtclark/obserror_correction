out = matrix(nrow = 1000, ncol =2)
for(i in 1:1000) {
  num = rnorm(2,4,1.2)
  out[i,1] = var(num)
  #sum((num-sum(num)/2)^2)
  #(2*(num[1]-num[2])^2)/4
  out[i,2] = ((num[1]-num[2])^2)/2
}
var(out)