%A=importdata('data.file',' ',0); %Funktioniert leider nicht bei der Größe!
%Import durch "Import Data"-Button! (als Matrix "A")
t=A(:,1);
alpha=A(:,2);
H=A(:,4);
alphadot=A(:,3);

figure(1) %Winkel in Abhängigkeit der Zeit
plot(t,alpha);

figure(2) %Energie
plot(t,H);

figure(3) %Winkel/Winkelgeschwindigkeit
plot(alpha,alphadot);

