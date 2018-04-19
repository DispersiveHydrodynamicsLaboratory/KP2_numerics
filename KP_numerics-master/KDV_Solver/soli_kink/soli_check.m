%% Figure out how to make a tanh that goes from ql to qr over a width w
ql = 1;
qr = 1/2;
w = 5;

a = (qr-ql)/2;
b = w;
c = 0;
d = a + ql;

x = -10:0.01:10;
y = (qr-ql)/2*tanh(w*(x-0)) + (qr+ql)/2;
y1 =             a*tanh(w*(x-c)) + d;

figure(1); 
    plot(x,y,x,y1);
