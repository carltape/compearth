function progreso(i,L)

if mod(i,ceil(L*0.01))==0
    p=round(i*100/L);
    disp([num2str(p) '% terminado']);
end
