p = 3;
for i=0:p
vector = load(['basis',num2str(i)]);
subplot(p+1,1,i+1);
plot(vector);
end
