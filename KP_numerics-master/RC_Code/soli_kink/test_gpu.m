soli = struct;
soli.th.intx = gpuArray(zeros(10,10));
soli.a = @(x) gpuArray(x);
soli.x = 1:10;
for Nxi =1:length(soli.x)
  s1=integral(@(x)gather(soli.a(x)),0,soli.x(Nxi),'ArrayValued',1);
  soli.th.intx(Nxi,:) = s1;
end
disp(soli.th.intx);
