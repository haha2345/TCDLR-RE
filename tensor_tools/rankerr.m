function result=rankerr(r,r_hat)
for i=1:3
result=abs(r_hat(i)-r(i))/3;
end