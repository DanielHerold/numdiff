T=1;
h=10^-5; %Zeiteinheit
y=zeros(2,round(T/h));
y(:,1)=[1;0]; %Anfangsposition der Katze
eps=10^-4;


for j=1:round(T/h)
    if norm([0;h*j]-y(:,j))<eps
        dead=j*h; %Gesamtschrittanzahl * Zeiteinheit
        break
    end
    y(:,j+1)=y(:,j)+h*2*([0;h*j]-y(:,j))/norm([0;h*j]-y(:,j));
end

disp( ['Die Katze erreicht die Maus nach ', num2str(dead),' Sekunden.' ]);
plot(y(1,:),y(2,:),'g');
axis([-0.1 1.1 -0.1 1.1]);
