b = 1;
y = exp(-(100/10000));
y1 = 2 * y * cos(2 * 3.142 * 890 * 0.0001)
ysq = y * y;
a = [1 y1 ysq]
for i = 1:300
  if((i == 1) | (i == 101) | (i == 201))
    x(i) = 1;
   else
    x(i) = 0;
   end
   if(i == 1)
    signal(i) = x(i);
   elseif(i == 2)
    signal(i) = x(i) + y1 * signal(i-1);
   else
    signal(i) = x(i) + y1 * signal(i-1) - ysq * signal(i-2);
   end
end
%sig = filter(b, a, x);
plot(signal);
