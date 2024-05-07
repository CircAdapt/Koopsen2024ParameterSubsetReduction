%% CircAdaptSetP
% Converts matlab P struct to C++ object
% At the moment, only single-patch Walmsley2015 version

function CA = CircAdaptSetP1(CA,P)
% disp('Set P to CA:');

%% General
% disp('-Set General');
nError = 0;
% nError = nError + 1 - CA.setScalar('','','','rhob',P.General.rhob);
nError = nError + 1 - CA.setScalar('','','','q0',P.General.q0);
nError = nError + 1 - CA.setScalar('','','','p0',P.General.p0);
%nError = nError + 1 - CA.setScalar('','','','tCycle',P.General.tCycle);
% nError = nError + 1 - CA.setScalar('','','','FacpControl',P.General.FacpControl);
% % something with ScaleVqY
% nError = nError + 1 - CA.setScalar('','','','dt',P.General.Dt);
% nError = nError + 1 - CA.setScalar('','','','tCycleRest',P.General.tCycleRest);
% nError = nError + 1 - CA.setScalar('','','','TimeFac',P.General.TimeFac);
% nError = nError + 1 - CA.setScalar('','','','PressFlowContr',P.General.PressFlowContr);
nError = nError + 1 - CA.setScalar('','','','dTauAv',P.General.dTauAv);

% if length(P.t)==1
%     nError = nError + 1 - CA.setScalar('','','','t',P.t);
% end

% disp(['-General parameters set with ' num2str(nError) ' error(s)']);

%% Solver
% if isfield(P,'Solver')
%     if isfield(P.Solver,'solverTol')
%         CA.setScalar('','','','solverTol',P.Solver.solverTol);
%     end
%     if isfield(P.Solver,'minTol')
%         CA.setScalar('','','','minTol',P.Solver.minTol);
%     end
%     if isfield(P.Solver,'runStableMaxBeats')
%         CA.setScalar('','','','maxBeats',P.Solver.runStableMaxBeats);
%     end
%     if isfield(P.Solver,'minmaxDt')
%         CA.setScalar('','','','minmaxDt',P.Solver.minmaxDt);
%     end
%     if isfield(P.Solver,'minDt')
%         CA.setScalar('','','','minDt',P.Solver.minDt);
%     end
%     if isfield(P.Solver,'maxLengthTimes1')
%         %CA.setScalar('','','','maxLengthTimes1',P.Solver.maxLengthTimes1);
%     end
% end

%% ArtVen
if isfield(P,'ArtVen') && ...
        P.ArtVen.n==2 
%     disp('-Set ArtVen');
    nError = 0;
    % Art and Ven properties
    %par = {'k','Len','A0','p0','AWall'};
    par = {'k','A0','p0'};
    for iP=1:length(par)
        nError = nError + 1 - CA.setScalar('CiSy','ArtVen','Art',par{iP},P.ArtVen.(par{iP})(1,1));
        nError = nError + 1 - CA.setScalar('CiSy','ArtVen','Ven',par{iP},P.ArtVen.(par{iP})(2,1));
%         nError = nError + 1 - CA.setScalar('CiPu','ArtVen','Art',par{iP},P.ArtVen.(par{iP})(1,2));
%         nError = nError + 1 - CA.setScalar('CiPu','ArtVen','Ven',par{iP},P.ArtVen.(par{iP})(2,2));
    end
    
    %artven properties
%     par = {'p0AV','q0AV','kAV'};
%     for iP=1:length(par)
%         nError = nError + 1 - CA.setScalar('CiSy','ArtVen','',par{iP},P.ArtVen.(par{iP})(1));
%         nError = nError + 1 - CA.setScalar('CiPu','ArtVen','',par{iP},P.ArtVen.(par{iP})(2));
%     end
%     disp(['-ArtVen parameters set with ' num2str(nError) ' error(s)']);
else
    error('ArtVen not compatible');
end

%% ArtVen
% if isfield(P,'Chamber') && ...
%         P.Chamber.n==2 
% %     disp('-All chamber parameter are stored in Cavity/Valve');
% else
%     error('ArtVen not compatible');
% end

%% TriSeg
if isfield(P,'TriSeg') && ...
        P.TriSeg.n==1
%     disp('-All TriSeg');
    nError = 0;
%     nError = nError + 1 - CA.setScalar('Lv','TriSeg','','tau',P.TriSeg.Tau(1));
    nError = nError + 1 - CA.setScalar('Lv','TriSeg','','V',P.TriSeg.V(end));
    nError = nError + 1 - CA.setScalar('Lv','TriSeg','','Y',P.TriSeg.Y(end));
%     disp(['-TriSeg parameters set with ' num2str(nError) ' error(s)']);
else
    error('TriSeg not compatible');
end

%% Cavity
if isfield(P,'Cavity') && ...
        P.Cavity.n==8
%     disp('-All Cavity');
    nError = 0;
    %artven properties
    loc1 = {'CiSy','CiSy','CiPu','CiPu','La','Ra','Lv','Rv'};
    loc2 = {'Art','Ven','Art','Ven','La','Ra','Lv','Rv'};
    for iL=1:length(loc1)
        nError = nError + 1 - CA.setScalar(loc1{iL},'Cavity',loc2{iL},'V',P.Cavity.V(end,iL));
    end
%     disp(['-Cavity parameters set with ' num2str(nError) ' error(s)']);
else
    error('Cavity not compatible');
end

%% Wall
% if isfield(P,'Wall') && ...
%         P.Wall.n==10
% %     disp('-All Wall');
%     nError = 0;
%     
%     loc1 = {'La','Ra','Lv','Sv','Rv'};
%     loc2 = {'','','','',''};
%     for iL=1:5
%         idx = strcmp(P.Wall.Name,loc1{iL});
%         nError = nError + 1 - CA.setScalar(loc1{iL},'Wall',loc2{iL},'AmDead',P.Wall.AmDead(end,idx));
%         
%         % check to split/merge patches
%         nPatch = CA.getScalar(loc1{iL},'Wall','','nPatch');
%         if P.Wall.nPatch(idx) < nPatch
%             % remove patches
% %             CA.setScalar(loc1{iL},'Wall',[loc1{iL} '1'],'mergePatch',nPatch - find( P.Wall.nPatch(idx) ) );
%         elseif P.Wall.nPatch(idx) > nPatch
%             % add patches
%             CA.setScalar(loc1{iL},'Wall',[loc1{iL} '1'],'Split',double(P.Wall.nPatch(idx) - nPatch + 1) );
%         end
%     end
% %     disp(['-Wall parameters set with ' num2str(nError) ' error(s)']);
% else
%     error('Wall not compatible');
% end

%% Patch
if isfield(P,'Patch')
%     disp('-All Patch');
    nError = 0;
    
    loc1 = {'La','Ra','Lv','Sv','Rv'};
    loc2 = {'La1','Ra1','Lv1','Sv1','Rv1'};
    
    loc1={};
    loc2={};
    for iP = 1:length(P.Patch.Name)
        loc1{iP} = P.Patch.Name{iP}(1:2);
        loc2{iP} = P.Patch.Name{iP};
    end
    
    %par = {'dT','Lsi','C','LsRef','Ls0Pas','dLsPas','SfPas','Lsi0Act','LenSeriesElement','SfAct','vMax','TimeAct','TR','TD','CRest','VWall','AmRef','ADO','LDAD','LDCI'};
    par = {'dT','Lsi','C','Ls0Pas','dLsPas','SfPas','LenSeriesElement','SfAct','vMax','AmRef','k1','VWall','TR','TD','ADO','LDAD','LDCI'};
        
%     if isfield(P.Patch,'k1')
%         par{end+1} = 'k1';
%     end
    
    for iL=1:length(loc1)
        for iP = 1:length(par)
            nError = nError + 1 - CA.setScalar(loc1{iL},'Patch',loc2{iL},par{iP},P.Patch.(par{iP})(end,iL));
        end
    end
%     disp(['-Patch parameters set with ' num2str(nError) ' error(s)']);
else
    error('Patch not compatible');
end

%% Node
%nothing

%% Bag

% if isfield(P,'Bag') && ...
%         P.Bag.n==1
% %     disp('-All Bag');
%     nError = 0;
%     nError = nError + 1 - CA.setScalar('Peri','Bag','','k',P.Bag.k(1));
%     nError = nError + 1 - CA.setScalar('Peri','Bag','','VRef',P.Bag.VRef(1));
%     nError = nError + 1 - CA.setScalar('Peri','Bag','','SfPeri',P.Bag.pAdapt(1));
% %     disp(['-Bag parameters set with ' num2str(nError) ' error(s)']);
% else
%     error('Bag not compatible');
% end

%% Valve
if isfield(P,'Valve') && ...
        P.Valve.n==9
%     disp('-All Valve, only first 6 valves');
    nError = 0;
    
    loc1 = P.Valve.Name;
    loc2 = loc1; % not needed here
    %par = {'q','AOpen','ALeak','Len'};
    par = {'q'};
    for iL=1:length(loc1)
        for iP = 1:length(par)
            nError = nError + 1 - CA.setScalar(loc1{iL},'Valve',loc2{iL},par{iP},P.Valve.(par{iP})(end,iL));
        end
    end
%     disp(['-Valve parameters set with ' num2str(nError) ' error(s)']);
else
    error('Valve not compatible');
end



%     
%     %% ArtVen
%     P.ArtVen.Name = {'Sy','Pu'};
%     P.ArtVen.n = 2;
%     P.ArtVen.iCavity = [1 3];
%     P.ArtVen.iWall = [1 3];
%     P.ArtVen.Adapt = struct;
%     %get Art and Ven data
%     CppLoc = {'CiSy','CiPu'};
%     CppLoc1 = {'Art','Ven'};
%     PPatchPar = {'k','Len','A0','p0','AWall'};
%     CppPar = PPatchPar;
%     PModule = 'ArtVen';
%     CppModule = 'ArtVen';
%     P=loadPmoduleScalarMatrix(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar);
%     %get ArtVen data
%     CppLoc1 = {'',''};
%     PPatchPar = {'p0AV','q0AV','kAV','q'};
%     CppPar = PPatchPar;
%     CppIsVec = [0 0 0 1];
%     P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
%        
%     %% Chamber
%     P.Chamber.Name = {'La','Ra'};
%     P.Chamber.n = 2;
%     P.Chamber.iCavity = [5 6];
%     P.Chamber.iWall = [5 6];
%     
%     %% TriSeg
%     P.TriSeg.Name = {'v'};
%     P.TriSeg.n = 1;
%     P.TriSeg.iCavity = 7;
%     P.TriSeg.iWall = 7;
%     P.TriSeg.V = CA.getVector('Lv','TriSeg','','V');
%     P.TriSeg.Y = CA.getVector('Lv','TriSeg','','Y');
%     P.TriSeg.Tau = CA.getScalar('Lv','TriSeg','','tau');
%     P.TriSeg.VS = CA.getVector('Lv','TriSeg','','VS');
%     P.TriSeg.YS = CA.getVector('Lv','TriSeg','','YS');
%     P.TriSeg.VDot = CA.getVector('Lv','TriSeg','','VDot');
%     P.TriSeg.YDot = CA.getVector('Lv','TriSeg','','YDot');
%     
%     
%     %% Cavity
%     P.Cavity.Name = {'SyArt','SyVen','PuArt','PuVen','La','Ra','Lv','Rv'};
%     P.Cavity.n = 8;
%     P.Cavity.iNode = [1 2 3 4 5 6 7 8];
%     P.Cavity.Adapt = struct;
%     % load data
%     CppLoc = {'CiSy','CiSy','CiPu','CiPu','La','Ra','Lv','Rv'};
%     CppLoc1 = {'Art','Ven','Art','Ven','','','',''};
%     PPatchPar = {'V','A','Z','p','VDot'};
%     CppPar = PPatchPar;
%     CppIsVec = [1 1 1 1 1];
%     PModule = 'Cavity';
%     CppModule = 'Cavity';
%     P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
%         
%     %% Wall
%     P.Wall.Name = {'SyArt','SyVen','PuArt','PuVen','La','Ra','Lv','Sv','Rv','Peri'};
%     P.Wall.n = 10;
%     P.Wall.nPatch = [0 0 0 0 1 1 1 1 1 0];
%     P.Wall.iPatch = [1 1 1 1 1 2 3 4 5 6];
%     %P.Wall.AmDead=
%     P.Wall.Am0 = [zeros(length(P.t),4) CA.getVector('La','Wall','','Am0') CA.getVector('Ra','Wall','','Am0') ...
%                                        CA.getVector('Lv','Wall','','Am0') CA.getVector('Sv','Wall','','Am0') CA.getVector('Rv','Wall','','Am0') zeros(length(P.t),1)];
%     P.Wall.DADT = [zeros(length(P.t),4) CA.getVector('La','Wall','','DADT') CA.getVector('Ra','Wall','','DADT') ...
%                                         CA.getVector('Lv','Wall','','DADT') CA.getVector('Sv','Wall','','DADT') CA.getVector('Rv','Wall','','DADT') zeros(length(P.t),1)];
%     P.Wall.T = [zeros(length(P.t),4) CA.getVector('La','Wall','','T') CA.getVector('Ra','Wall','','T') ...
%                                      CA.getVector('Lv','Wall','','T') CA.getVector('Sv','Wall','','T') CA.getVector('Rv','Wall','','T') zeros(length(P.t),1)];
%     P.Wall.Cm = [zeros(length(P.t),4) CA.getVector('La','Wall','','Cm') CA.getVector('Ra','Wall','','Cm') ...
%                                       CA.getVector('Lv','Wall','','Cm') CA.getVector('Sv','Wall','','Cm') CA.getVector('Rv','Wall','','Cm') zeros(length(P.t),1)];
%     P.Wall.Am = [zeros(length(P.t),4) CA.getVector('La','Wall','','Am') CA.getVector('Ra','Wall','','Am') ...
%                                       CA.getVector('Lv','Wall','','Am') CA.getVector('Sv','Wall','','Am') CA.getVector('Rv','Wall','','Am') zeros(length(P.t),1)];
%     P.Wall.pTrans = [zeros(length(P.t),4) CA.getVector('La','Wall','','pTrans') CA.getVector('Ra','Wall','','pTrans') ...
%                                           CA.getVector('Lv','Wall','','pTrans') CA.getVector('Sv','Wall','','pTrans') CA.getVector('Rv','Wall','','pTrans') zeros(length(P.t),1)];
%     % P.Wall.VWall = 
%     
%     %% Patch
%     P.Patch.Name = {'La1','Ra1','Lv1','Sv1','Rv1'};
%     P.Patch.n = 5;
%     
%     % load patch before adapt
%     CppLoc = {'La','Ra','Lv','Sv','Rv'};
%     CppLoc1 = P.Patch.Name;
%     PPatchPar = {'dT','Lsi','C','ActivationDelay','LsRef','Ls0Pas','dLsPas','SfPas','Lsi0Act','LenSeriesElement','SfAct','vMax','TimeAct','TR','TD','CRest','VWall','AmRef'};
%     CppPar = PPatchPar;
%     CppIsVec = [0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%     PModule = 'Patch';
%     CppModule = 'Patch';
%     P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
%    
%     P.Patch.Adapt = struct;
%     
%     %load patch after adapt
%     CppLoc = {'La','Ra','Lv','Sv','Rv'};
%     CppLoc1 = P.Patch.Name;
%     PPatchPar = {'T','Ef','Ls','CDot','SfEcm','SfPasT','LsiDot','Sf','DSfDEf','DADT','Am0','Am'};
%     CppPar = PPatchPar;
%     CppIsVec = [1 1 1 1 1 1 1 1 1 1 1 1 1];
%     PModule = 'Patch';
%     CppModule = 'Patch';
%     P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
%   
%      % Node
%     P.Node.Name = {'SyArt','SyVen','PuArt','PuVen','La','Ra','Lv','Rv'};
%     P.Node.n = 8;
%     P.Node.iCavity = [1 2 3 4 5 6 7 8];
%     P.Node.q = [CA.getVector('CiSy','ArtVen','Art','q') CA.getVector('CiSy','ArtVen','Ven','q') ...
%                   CA.getVector('CiPu','ArtVen','Art','q') CA.getVector('CiPu','ArtVen','Ven','q') ...
%                   CA.getVector('La','Chamber','La1','q') CA.getVector('Ra','Chamber','Ra1','q') ...
%                   CA.getVector('Lv','Cavity','','q') CA.getVector('Rv','Cavity','','q')];
%     P.Node.Y = [CA.getVector('CiSy','ArtVen','Art','Y') CA.getVector('CiSy','ArtVen','Ven','Y') ...
%                   CA.getVector('CiPu','ArtVen','Art','Y') CA.getVector('CiPu','ArtVen','Ven','Y') ...
%                   CA.getVector('La','Chamber','La1','Y') CA.getVector('Ra','Chamber','Ra1','Y') ...
%                   CA.getVector('Lv','Cavity','','Y') CA.getVector('Rv','Cavity','','Y')];
%     P.Node.p = [CA.getVector('CiSy','ArtVen','Art','p') CA.getVector('CiSy','ArtVen','Ven','p') ...
%                   CA.getVector('CiPu','ArtVen','Art','p') CA.getVector('CiPu','ArtVen','Ven','p') ...
%                   CA.getVector('La','Chamber','La1','p') CA.getVector('Ra','Chamber','Ra1','p') ...
%                   CA.getVector('Lv','Cavity','','p') CA.getVector('Rv','Cavity','','p')];
%               
%     % Bag
%     P.Bag.Name = {'Peri'};
%     P.Bag.n = 1;
%     P.Bag.iWall = 10;
% %     P.Bag.VRef;
%     P.Bag.p = CA.getVector('Peri','Bag','','p');
%     
%     %% Valve
%     P.Valve.Name = {'SyVenRa','RaRv','RvPuArt','PuVenLa','LaLv','LvSyArt','LaRa','LvRv','SyArtPuArt'};
%     P.Valve.n = 9;
%     P.Valve.iNodeProx = [2 6 8 4 5 7 5 7 1];
%     P.Valve.iNodeDist = [6 8 3 5 7 1 6 8 3];
%     
%     P.Valve.Name = {'SyVenRa','RaRv','RvPuArt','PuVenLa','LaLv','LvSyArt'};
%     P.Valve.n = 6;
%     P.Valve.iNodeProx = [2 6 8 4 5 7];
%     P.Valve.iNodeDist = [6 8 3 5 7 1];
%     
%     CppLoc = {P.Valve.Name{1:6}} ;
%     CppLoc1 = {'','','','','',''};
%     PPatchPar = {'q','AOpen','ALeak','Len','L','A','qDot'};
%     CppPar = PPatchPar;
%     CppIsVec = [1 0 0 0 1 1 1];
%     PModule = 'Valve';
%     CppModule = 'Valve';
%     P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
end

function P = loadPmoduledata(P,CA,Pmodule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec)
    for iP = 1:length(CppPar)
        if CppIsVec(iP)
            curPar = zeros(length(P.t),length(CppLoc));
            for iL=1:length(CppLoc)
                curPar(:,iL) = CA.getVector(CppLoc{iL},CppModule,CppLoc1{iL},CppPar{iP});
            end
        else
            curPar = zeros(1,length(CppLoc));
            for iL=1:length(CppLoc)
                curPar(iL) = CA.getScalar(CppLoc{iL},CppModule,CppLoc1{iL},CppPar{iP});
            end
        end
        P.(Pmodule).(PPatchPar{iP}) = curPar;
    end
end

% Designed for ArtVen
function P = loadPmoduleScalarMatrix(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar)
    for iP = 1:length(CppPar)
        curPar = zeros(length(CppLoc1),length(CppLoc));
        for iL=1:length(CppLoc)
            for iL1=1:length(CppLoc1)
                curPar(iL1,iL) = CA.getScalar(CppLoc{iL},CppModule,CppLoc1{iL1},CppPar{iP});
            end
        end
        P.(PModule).(PPatchPar{iP}) = curPar;
    end
end