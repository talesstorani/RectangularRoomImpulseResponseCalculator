function [RoomData]=ImageSourceRectangularRoomImpulseResponseCalculator(samplingFrequency,receiverPosition ,sourcePosition , roomDimensions ,reflexionOrder, frequencyIntervalDesired ,absorptionCoefficients,temperature,signalToNoiseRatio)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Function for calculation of Impulse response, frequency response     %
%    function and average reverberation time of rectangular rooms         %
%    using the image source model.                                        %
%                                                                         %
%   Author: T. MARINI STORANI                                             %
%   Version 1.1 - 20/05/2018                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% output
% RoomData (struct) entrega 3 parametros
% RoomData.ImpulseResponse = resposta impulsiva 'symmetric'
% RoomData.FrequencyResponse = Fun?ao de resposta em frequ?ncia
% RoomData.ReverberationTime = tempo de reverbera??o m?dio calculado pela rotina 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         COORDINATES                ^                    %   
%                                                    |                    %   
%                                                    |                    %     
%                                                    |y2                  %        
%            _ _ _                                   |                    %                                                      
%          |\     \                   ^     _ _ _ _ _|                    %        
%          | \  3  \                  |    /    x2   /                    %  
%          |  \_ _ _\              y1 |   /         /                     %
%          \ 2|     |                 |  /z1       /z2                    %
%        z  \ |  1  | y               | /         /                       % 
%            \o_ _ _|                 |/         /                        %         
%                x                    o_ _ _ _ _/                         %
%             Fig (1)                      x1                             %
%                                        Fig (2)                          %
%                                                                         %
%                   origem  ->  o=(0 , 0, 0)                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%                                 Read me.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algoritmo ? pensado na forma dos artigos: 
% "Image method for efficiently simulating small-room acoustics" de 
% Jont B. Allen and David A. Berkley
%
% 
% Uses function for permutation with repetition from:
% "Jos (10584) (2022). permn (https://www.mathworks.com/matlabcentral/fileexchange/7147-permn), MATLAB Central File Exchange. Retrieved August 29, 2022."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AS COORDENADAS COM INDICE "1" S?O REFERENTES AS PAREDES QUE CORTAM O 
% PONTO "o" DA ORIGEM. AS COORDENADAS DO ?NDICE "2" S?O REFERENTES AS
% PAREDES QUE N?O PASSAM PELA ORIGEM. AS ?NDICES S?O USADOS PARA
% LOCALIZA??O E COLOCA??O DA ABSOR??O CORRETA NO VETOR ABSOR??O. POR EXEMP. 
% alpha = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6] ENT?O ESTE VETOR SE REFERE
% alpha = [x1 , x2 , y1 , y2 , z1 , z2 ] AS PERMUTA??ES q,j,k, ENTRE 0-1 e 
% n,l,m, ENTRE -N a +N, INDICAR?O QUAIS PAREDES O RAIO EMITIDO PELAS
% FONTE-IMAGENS ATINGIR?O.  N?MEROS DAS PAREDES S?O IGUAIS AOS 
% DADOS(de jogos de tabeileiro),OU SEJA, SOMA ENTRE LADOS PARALELOS DEVE 
% SER IGUAL A 7, vide Fig. (1). ENT?O, OS COEFICIENTES DE ABSOR??O alpha  
% DEVEM SER COLOCADOS CONFORME A ORDEM DE PAREDES DADA POR: 
% alpha=[1,3,6,5,2,4], EM QUE OS NUMEROS DENTRO DO VETOR alpha CORRESPONDEM 
% AOS COEFICIENTES DE ABSOR??O M?DIO DE CADA PAREDE SEGUINDO ORDEM Fig(1).


% Rotina RI_sala_retang calcula resposta impulsiva de salas retangulares 
% em banda de oitava, ter?o de oitava e estreita. Quanto maior o tempo de 
% reverbera??o, maior ordem de reflex?o(n) e maior n? bandas de frequ?ncia
% maior ser? o tempo de processamento. 


%% EXEMPLE DE PAR?METROS INICIAIS 

% samplingFrequency=24000;                
% % % [ x  ,  y  , z]
% sourcePosition=  [3.7 , 1.92 , 1.3];    
% receiverPosition=[5.89   , 1.6  , 3.5];    
% roomDimensions=[7.37 , 3 , 9.85];    
% reflexionOrder=10;                    
% frequencyIntervalDesired=[20 11000];         
% absorptionCoefficients=[0.06 , 0.06 , 0.06 , 0.06 , 0.06 , 0.27];  
% temperature=23;                
% signalToNoiseRatio=60;           

%% Room Parameters
 
RectRoomParameters.SourcePosition=sourcePosition;     
RectRoomParameters.ReceiverPosition=receiverPosition;      
RectRoomParameters.RoomDimensions=roomDimensions;       
RectRoomParameters.ReflexionOrder=reflexionOrder;                
RectRoomParameters.AbsorptionCoefficients=absorptionCoefficients;
RectRoomParameters.Temperature=temperature;                   
RectRoomParameters.SignalToNoiseRatio=signalToNoiseRatio;                 

%% TESTES ANTI ERRO

% TESTA PARA QUANTIDADE DE ALPHA/totaldePAREDEs INCOERENTE (RAZ?O DIFERENTE DE 1)

if  length(absorptionCoefficients)/6 ~= 1
    error('Quantidade de elementos de absor?ao est? maior ou menor que o n?mero total de paredes');
end

% TESTA PARA VALORES INCOERENTES DE DIMENS?ES, LOCALIZA??O E ABSOR??O
if RectRoomParameters.SourcePosition(1)>=RectRoomParameters.RoomDimensions(1) || RectRoomParameters.SourcePosition(2)>=RectRoomParameters.RoomDimensions(1) || RectRoomParameters.SourcePosition(1)>=RectRoomParameters.RoomDimensions(3) || RectRoomParameters.ReceiverPosition(1)>=RectRoomParameters.RoomDimensions(1)...
   || RectRoomParameters.ReceiverPosition(2)>=RectRoomParameters.RoomDimensions(2) || RectRoomParameters.ReceiverPosition(3)>=RectRoomParameters.RoomDimensions(3) 
error('Fonte ou Microfone est?o fora das dimens?es da sala');
end

% TESTA PARA VALORES DE ALPHA INCOERENTES;
if any(find(RectRoomParameters.AbsorptionCoefficients(:)>1)) || any(find(RectRoomParameters.AbsorptionCoefficients(:)<0))
        error('Absor??o deve estar contida no intervalo 0 - 1');
end

% TESTA PARA FREQUENCIA AMOSTRAL
if samplingFrequency<2*frequencyIntervalDesired(end)
error('frequ?ncia amostral n?o pode ser menor que duas vezes a m?xima freq?ncia de an?lise');
end

reflexionCoefficients=abs((sqrt((RectRoomParameters.AbsorptionCoefficients)-1)));

%% CALCULOS DOS TR's e Tmax para avalia??o da RI

% % VELOCIDADE DO SOM
RectRoomParameters.soundSpeed=20.05*sqrt(273+RectRoomParameters.Temperature);

% T60 m?dio EYRING DA SALA   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dimensionsForReverberationTime= permn(RectRoomParameters.RoomDimensions,2); 
iteratorOverDimension=1:4:length(dimensionsForReverberationTime);
if dimensionsForReverberationTime(iteratorOverDimension,1)==dimensionsForReverberationTime(iteratorOverDimension,2)
  dimensionsForReverberationTime(iteratorOverDimension,:)=[];   
end
dimIndex=1:6;
dimensionsForReverberationTime=dimensionsForReverberationTime(dimIndex,1).*dimensionsForReverberationTime(dimIndex,2);

RectRoomParameters.ReverberationTime=((4.*log(10e-6).*prod(RectRoomParameters.RoomDimensions)))./(RectRoomParameters.soundSpeed.*(sum(dimensionsForReverberationTime)).*log(1-((sum(dimensionsForReverberationTime'.*RectRoomParameters.AbsorptionCoefficients)./(sum(dimensionsForReverberationTime))))));
% pode-se usar o de Sabine. Descomenta-se embaixo e comenta em cima
% RectRoom.ReverberationTime(i)=(24*log(10)*prod(RectRoom.RoomDimensions))/(RectRoom.soundSpeed*sum(DIMtr'.*absor(i).s));
RectRoomParameters.TemperatureAndSNR=RectRoomParameters.SignalToNoiseRatio.*RectRoomParameters.ReverberationTime./60;
RectRoomParameters.MaxTimeLength=RectRoomParameters.TemperatureAndSNR+(norm(RectRoomParameters.SourcePosition-RectRoomParameters.ReceiverPosition)/RectRoomParameters.soundSpeed);
RectRoomParameters.MaxOrder = (ceil(RectRoomParameters.MaxTimeLength*samplingFrequency));

%%  PERMUTA??ES DAS DIMENS?ES DA SALA
% realiza todas as permuta??es dos valores de ordem de reflex?o
% resultado ? uma matriz (2N+1)^3
% perm=struct de todas as permuta??es realizadas no c?digo
perm.NthPermutation=(permn(-RectRoomParameters.ReflexionOrder:RectRoomParameters.ReflexionOrder,3));

%% MATRIZ PERMUTA??ES DO TOTAL DE PAREDES E TAMANHO DA SALA  
perm.RoomWallPermutation=(2.*permn(-RectRoomParameters.ReflexionOrder:RectRoomParameters.ReflexionOrder,3)).*repmat(RectRoomParameters.RoomDimensions,length(perm.NthPermutation),1);

%% MATRIZ PERMUTA??ES DOS SUB ?NDICES q j k PARA USO EM MIC1
perm.subIndexesPermutations=permn([0 1],3);

%% MATRIZ REPETI??O DO VETOR MIC PARA ACERTAR A MATRIX DE PERMUTA??O 
RectRoomParameters.ReceiverPosition1=((2*perm.subIndexesPermutations)-1).*repmat(RectRoomParameters.ReceiverPosition,8,1);

%% RAIO SOMA F-MIC PARA TODAS AS PERMUTA??ES DENTRO DA SALA E CADA PAREDE
perm.SumSourceReceiverPermutation=repmat(RectRoomParameters.SourcePosition,8,1)+RectRoomParameters.ReceiverPosition1;

%% PERMUTA??O DE CADA VALOR INTEIRO DE SUB ?NDICE DA MATRIZ DE RAIO 
% PARA OBTEN??O DA MATRIZ DE TODOS OS RAIOS DE DIST?NCIA ENTRE (F -MIC)
perm.RoomWallPermutation=[repmat(perm.SumSourceReceiverPermutation(1,:),length(perm.RoomWallPermutation),1)+ perm.RoomWallPermutation(1:length(perm.RoomWallPermutation),:)...
        ;repmat(perm.SumSourceReceiverPermutation(2,:),length(perm.RoomWallPermutation),1)+ perm.RoomWallPermutation(1:length(perm.RoomWallPermutation),:)...
        ;repmat(perm.SumSourceReceiverPermutation(3,:),length(perm.RoomWallPermutation),1)+ perm.RoomWallPermutation(1:length(perm.RoomWallPermutation),:)...
        ;repmat(perm.SumSourceReceiverPermutation(4,:),length(perm.RoomWallPermutation),1)+ perm.RoomWallPermutation(1:length(perm.RoomWallPermutation),:)...
        ;repmat(perm.SumSourceReceiverPermutation(5,:),length(perm.RoomWallPermutation),1)+ perm.RoomWallPermutation(1:length(perm.RoomWallPermutation),:)...
        ;repmat(perm.SumSourceReceiverPermutation(6,:),length(perm.RoomWallPermutation),1)+ perm.RoomWallPermutation(1:length(perm.RoomWallPermutation),:)...
        ;repmat(perm.SumSourceReceiverPermutation(7,:),length(perm.RoomWallPermutation),1)+ perm.RoomWallPermutation(1:length(perm.RoomWallPermutation),:)...
        ;repmat(perm.SumSourceReceiverPermutation(8,:),length(perm.RoomWallPermutation),1)+ perm.RoomWallPermutation(1:length(perm.RoomWallPermutation),:)];

%% EXPOENTE DAS REFLEX?ES EM FUN??O DOS SUB ?NDICES DA MATRIZ DE RAIO
perm.nq=repmat(perm.NthPermutation(:,1),length(perm.subIndexesPermutations),1)+repmat(perm.subIndexesPermutations(:,1),length(perm.NthPermutation),1);
perm.n=repmat(perm.NthPermutation(:,1),length(perm.subIndexesPermutations),1);
perm.lj=repmat(perm.NthPermutation(:,2),length(perm.subIndexesPermutations),1)+repmat(perm.subIndexesPermutations(:,2),length(perm.NthPermutation),1);
perm.l=repmat(perm.NthPermutation(:,2),length(perm.subIndexesPermutations),1);
perm.mk=repmat(perm.NthPermutation(:,3),length(perm.subIndexesPermutations),1)+repmat(perm.subIndexesPermutations(:,3),length(perm.NthPermutation),1);
perm.m=repmat(perm.NthPermutation(:,3),length(perm.subIndexesPermutations),1);

%% Radius of every reflexion order
%GEOMETRIA EUCLIDIANA PARA CALCULO DO VETOR DE DISTANCIAS ENTRE F e MIC
% obs: MATLAB 2017+ possui fun??o espec?fica e pode ser substituido
radius=sqrt(sum(perm.RoomWallPermutation.^2, 2));

%% CALCULO DOS PASSOS DE FREQUENCIAS
% passos de frequencias est?o em potencia de 2. Para maior rapidez.
%testa para ordem de reflex?o dos atrasos menores que o TR. Usu?rio pode
%seguir com o c?digo, mas n?o ser? t?o representativo.
% No passo de frequencia se o atraso for maior que o TR usa ordem maxima
% para claculo pois ? necess?rio um vetor de amostras maior para melhor
% representatividade. Sen?o ele usa o espa?o entre freq(i)inicial 
% at? freq(i) final tamb?m espa?ado em um multiplo de pot?ncia de 2. 
    
  %PASSO DE FREQUENCIA BANDA ESTREITA  
if RectRoomParameters.ReverberationTime< max(RectRoomParameters.ReverberationTime)/2.5        
   frequencia=2*pi*linspace(frequencyIntervalDesired(1),frequencyIntervalDesired(2),2^nextpow2(length(frequencyIntervalDesired(1):frequencyIntervalDesired(2))));
else
   frequencia=2*pi*linspace(frequencyIntervalDesired(1),frequencyIntervalDesired(2),2^nextpow2(RectRoomParameters.MaxOrder));
end

%% COEFICIENTE DE REFLEX?O EM FUN??O DA ABSOR??O ALPHA

reflexionCoefficientsArray=reflexionCoefficients(1).^(abs(perm.nq)) .* reflexionCoefficients(2).^(abs(perm.n))...
                        .* reflexionCoefficients(3).^(abs(perm.lj)) .* reflexionCoefficients(4).^(abs(perm.l))...  
                        .* reflexionCoefficients(5).^(abs(perm.mk)) .* reflexionCoefficients(6).^(abs(perm.m));

%% CALCULA TODAS AS INTERA??ES ENTRE FREQ E RAIO PARA CALCULO DA FRF Hw

tic

frequencyResponse=zeros(1,length(frequencia));
for rediusIterator=1:length(radius)    
    totalComplexSoundPressureEnergy=(reflexionCoefficientsArray(rediusIterator)./(4.*pi.*radius(rediusIterator))).*exp(-1i.*frequencia.*radius(rediusIterator)./RectRoomParameters.soundSpeed);%      
    frequencyResponse=totalComplexSoundPressureEnergy+frequencyResponse;    
end

impulseResponse=real(ifft(frequencyResponse));
impulseResponse=impulseResponse(1:length(frequencyResponse)/2);

%% Final struct

RoomData.ImpulseResponse=impulseResponse;
RoomData.FrequencyResponse=frequencyResponse(1:length(frequencyResponse));
RoomData.ReverberationTime.Value=RectRoomParameters.ReverberationTime;
RoomData.ReverberationTime.CalculationType='Eyring'; 
RoomData.InitialInputData=RectRoomParameters;
RoomData.TotalCalculationTime=toc;

end

                                                               