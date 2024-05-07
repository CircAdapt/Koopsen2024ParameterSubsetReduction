%% CircAdaptGetP
% Converts C++ object to matlab P struct

function P = CircAdaptGetP(CA,P)
    % init
    if iscell(P) && length(P)==1
        P = P{1};
    elseif length(P)>1
        error('niet toegestaan');
    end
    
    if ~exist('P','var') || isempty(P)
        P = struct;
        P.General = struct;
        P.ArtVen = struct;
        P.Chamber = struct;
        P.TriSeg = struct;
        P.Cavity = struct;
        P.Wall = struct;
        P.Patch = struct;
        P.Node = struct;
        P.Bag = struct;
        P.Valve = struct;
        P.Adapt = struct;
    end
    
    P.t = CA.getVector('','','','Time');
    
    %% Solver Data
    P.Log.Code    = 'cpp';
    P.Log.Version = CA.getVersion();
    
    %% Solver
    P.Solver.solverTol         = CA.getScalar('','','','solverTol');
    P.Solver.minTol            = CA.getScalar('','','','minTol');
    P.Solver.runStableMaxBeats = CA.getScalar('','','','maxBeats');
    P.Solver.minmaxDt          = CA.getScalar('','','','minmaxDt');
    P.Solver.minDt             = CA.getScalar('','','','minDt');
    P.Solver.maxLengthTimes1   = CA.getScalar('','','','maxLengthTimes1');

    
    %% General
    P.General.rhob = CA.getScalar('','','','rhob');
    P.General.q0 = CA.getScalar('','','','q0');
    P.General.p0 = CA.getScalar('','','','p0');
    P.General.tCycle = CA.getScalar('','','','tCycle');
    P.General.FacpControl = CA.getScalar('','','','FacpControl');
    P.General.ScaleVqY = [1e-05,0.0001,0.1];
    P.General.Dt = CA.getScalar('','','','dt');
    P.General.tCycleRest = CA.getScalar('','','','tCycleRest');
    P.General.TimeFac = CA.getScalar('','','','TimeFac');
    P.General.PressFlowContr = CA.getScalar('','','','PressFlowContr');
    P.General.TauAv = CA.getScalar('','','','TauAv');
    P.General.dTauAv = CA.getScalar('','','','dTauAv');

    
    %% ArtVen
    P.ArtVen.Name = {'Sy','Pu'};
    P.ArtVen.n = 2;
    P.ArtVen.iCavity = [1 3];
    P.ArtVen.iWall = [1 3];
    P.ArtVen.Adapt = struct;
    %get Art and Ven data
    CppLoc = {'CiSy','CiPu'};
    CppLoc1 = {'Art','Ven'};
    PPatchPar = {'k','Len','A0','p0','AWall'};
    CppPar = PPatchPar;
    PModule = 'ArtVen';
    CppModule = 'ArtVen';
    P=loadPmoduleScalarMatrix(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar);
    %get ArtVen data
    CppLoc1 = {'',''};
    PPatchPar = {'p0AV','q0AV','kAV','q'};
    CppPar = PPatchPar;
    CppIsVec = [0 0 0 1];
    P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
       
    %% Chamber
    P.Chamber.Name = {'La','Ra'};
    P.Chamber.n = 2;
    P.Chamber.iCavity = [5 6];
    P.Chamber.iWall = [5 6];
    
    %% TriSeg
    P.TriSeg.Name = {'v'};
    P.TriSeg.n = 1;
    P.TriSeg.iCavity = 7;
    P.TriSeg.iWall = 7;
    P.TriSeg.V = CA.getVector('Lv','TriSeg','','V');
    P.TriSeg.Y = CA.getVector('Lv','TriSeg','','Y');
    P.TriSeg.Tau = CA.getScalar('Lv','TriSeg','','tau');
    P.TriSeg.VS = CA.getVector('Lv','TriSeg','','VS');
    P.TriSeg.YS = CA.getVector('Lv','TriSeg','','YS');
    P.TriSeg.VDot = CA.getVector('Lv','TriSeg','','VDot');
    P.TriSeg.YDot = CA.getVector('Lv','TriSeg','','YDot');
    
    
    %% Cavity
    P.Cavity.Name = {'SyArt','SyVen','PuArt','PuVen','La','Ra','Lv','Rv'};
    P.Cavity.n = 8;
    P.Cavity.iNode = [1 2 3 4 5 6 7 8];
    P.Cavity.Adapt = struct;
    % load data
    CppLoc = {'CiSy','CiSy','CiPu','CiPu','La','Ra','Lv','Rv'};
    CppLoc1 = {'Art','Ven','Art','Ven','','','',''};
    PPatchPar = {'V','A','Z','p','VDot'};
    CppPar = PPatchPar;
    CppPar{4} = 'pCavity';
    CppIsVec = [1 1 1 1 1];
    PModule = 'Cavity';
    CppModule = 'Cavity';
    P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
        
    %% Wall
    P.Wall.Name = {'SyArt','SyVen','PuArt','PuVen','La','Ra','Lv','Sv','Rv','Peri'};
    P.Wall.n = 10;
    P.Wall.nPatch = [0 0 0 0 1 1 1 1 1 0];
    P.Wall.iPatch = [1 1 1 1 1 2 3 4 5 6];
    %P.Wall.AmDead=
    P.Wall.AmDead = [zeros(1,4) CA.getScalar('La','Wall','','AmDead') CA.getScalar('Ra','Wall','','AmDead') ...
                                       CA.getScalar('Lv','Wall','','AmDead') CA.getScalar('Sv','Wall','','AmDead') CA.getScalar('Rv','Wall','','AmDead') zeros(1,1)];
    P.Wall.VWall = [zeros(1,4) CA.getScalar('La','Wall','','VWall') CA.getScalar('Ra','Wall','','VWall') ...
                                       CA.getScalar('Lv','Wall','','VWall') CA.getScalar('Sv','Wall','','VWall') CA.getScalar('Rv','Wall','','VWall') zeros(1,1)];
    P.Wall.Am0 = [zeros(length(P.t),4) CA.getVector('La','Wall','','Am0') CA.getVector('Ra','Wall','','Am0') ...
                                       CA.getVector('Lv','Wall','','Am0') CA.getVector('Sv','Wall','','Am0') CA.getVector('Rv','Wall','','Am0') zeros(length(P.t),1)];
    P.Wall.DADT = [zeros(length(P.t),4) CA.getVector('La','Wall','','DADT') CA.getVector('Ra','Wall','','DADT') ...
                                        CA.getVector('Lv','Wall','','DADT') CA.getVector('Sv','Wall','','DADT') CA.getVector('Rv','Wall','','DADT') zeros(length(P.t),1)];
    P.Wall.T = [zeros(length(P.t),4) CA.getVector('La','Wall','','T') CA.getVector('Ra','Wall','','T') ...
                                     CA.getVector('Lv','Wall','','T') CA.getVector('Sv','Wall','','T') CA.getVector('Rv','Wall','','T') zeros(length(P.t),1)];
    P.Wall.Cm = [zeros(length(P.t),4) CA.getVector('La','Wall','','Cm') CA.getVector('Ra','Wall','','Cm') ...
                                      CA.getVector('Lv','Wall','','Cm') CA.getVector('Sv','Wall','','Cm') CA.getVector('Rv','Wall','','Cm') zeros(length(P.t),1)];
    P.Wall.Am = [zeros(length(P.t),4) CA.getVector('La','Wall','','Am') CA.getVector('Ra','Wall','','Am') ...
                                      CA.getVector('Lv','Wall','','Am') CA.getVector('Sv','Wall','','Am') CA.getVector('Rv','Wall','','Am') zeros(length(P.t),1)];
    P.Wall.Vm = [zeros(length(P.t),4) CA.getVector('La','Wall','','Vm') CA.getVector('Ra','Wall','','Vm') ...
                                      CA.getVector('Lv','Wall','','Vm') CA.getVector('Sv','Wall','','Vm') CA.getVector('Rv','Wall','','Vm') zeros(length(P.t),1)];
    P.Wall.pTrans = [zeros(length(P.t),4) CA.getVector('La','Wall','','pTrans') CA.getVector('Ra','Wall','','pTrans') ...
                                          CA.getVector('Lv','Wall','','pTrans') CA.getVector('Sv','Wall','','pTrans') CA.getVector('Rv','Wall','','pTrans') zeros(length(P.t),1)];
    % P.Wall.VWall = 
    
    %% Patch
%     P.Patch.Name = {'La1','Ra1','Lv1','Sv1','Rv1'};
%     P.Patch.n = 5;
%     
    P.Patch.Name = {};
%     
    wallnames = {'La','Ra','Lv','Sv','Rv'};
    CppLoc = {};
    for iW = 1:length(wallnames)
        nPatch = CA.getScalar(wallnames{iW},'Wall','','nPatch');
        P.Wall.nPatch(4+iW) = nPatch;
        for iP = 1:nPatch
            P.Patch.Name{end+1} = [wallnames{iW} num2str(iP)];
            CppLoc{end+1} = [wallnames{iW}];
        end
    end
    
    P.Patch.n = length(P.Patch.Name);
    
    % load patch before adapt
    CppLoc1 = P.Patch.Name;
    PPatchPar = {'dT','Lsi','C','ActivationDelay','LsRef','Ls0Pas','dLsPas','SfPas','k1','Lsi0Act','LenSeriesElement','SfAct','vMax','TimeAct','TR','TD','CRest','VWall','AmRef','ADO','LDAD','LDCI'};
    CppPar = PPatchPar;
    CppIsVec = [0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    PModule = 'Patch';
    CppModule = 'Patch';
    P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
    P.Patch.ActivationDelay = P.Patch.ActivationDelay + [-P.General.tCycle 0]';
    
    P.Patch.Adapt = struct;
    
    %load patch after adapt
    PPatchPar = {'T','Ef','Ls','CDot','SfEcm','SfPasT','LsiDot','Sf','DSfDEf','DADT','Am0','Am'};
    CppPar = PPatchPar;
    CppIsVec = [1 1 1 1 1 1 1 1 1 1 1 1 1];
    PModule = 'Patch';
    CppModule = 'Patch';
    P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
  
     % Node
    P.Node.Name = {'SyArt','SyVen','PuArt','PuVen','La','Ra','Lv','Rv'};
    P.Node.n = 8;
    P.Node.iCavity = [1 2 3 4 5 6 7 8];
    P.Node.q = [CA.getVector('CiSy','ArtVen','Art','q') CA.getVector('CiSy','ArtVen','Ven','q') ...
                  CA.getVector('CiPu','ArtVen','Art','q') CA.getVector('CiPu','ArtVen','Ven','q') ...
                  CA.getVector('La','Chamber','La1','q') CA.getVector('Ra','Chamber','Ra1','q') ...
                  CA.getVector('Lv','Cavity','','q') CA.getVector('Rv','Cavity','','q')];
    P.Node.Y = [CA.getVector('CiSy','ArtVen','Art','Y') CA.getVector('CiSy','ArtVen','Ven','Y') ...
                  CA.getVector('CiPu','ArtVen','Art','Y') CA.getVector('CiPu','ArtVen','Ven','Y') ...
                  CA.getVector('La','Chamber','La1','Y') CA.getVector('Ra','Chamber','Ra1','Y') ...
                  CA.getVector('Lv','Cavity','','Y') CA.getVector('Rv','Cavity','','Y')];
    P.Node.p = [CA.getVector('CiSy','ArtVen','Art','p') CA.getVector('CiSy','ArtVen','Ven','p') ...
                  CA.getVector('CiPu','ArtVen','Art','p') CA.getVector('CiPu','ArtVen','Ven','p') ...
                  CA.getVector('La','Chamber','La1','p') CA.getVector('Ra','Chamber','Ra1','p') ...
                  CA.getVector('Lv','Cavity','','p') CA.getVector('Rv','Cavity','','p')];
              
    % Bag
    P.Bag.Name = {'Peri'};
    P.Bag.n = 1;
    P.Bag.iWall = 10;
    P.Bag.pAdapt = CA.getScalar('Peri','Bag','','SfPeri');
    P.Bag.p = CA.getVector('Peri','Bag','','p');
    P.Bag.k = CA.getScalar('Peri','Bag','','k');
    P.Bag.VRef = CA.getScalar('Peri','Bag','','VRef');
    
    %% Valve
    P.Valve.Name = {'SyVenRa','RaRv','RvPuArt','PuVenLa','LaLv','LvSyArt','LaRa','LvRv','SyArtPuArt'};
    P.Valve.n = 9;
    P.Valve.iNodeProx = [2 6 8 4 5 7 5 7 1];
    P.Valve.iNodeDist = [6 8 3 5 7 1 6 8 3];
    
    
    CppLoc = P.Valve.Name ;
    CppLoc1 = {'','','','','','','','',''};
    PPatchPar = {'q','AOpen','ALeak','Len','L','A','qDot'};
    CppPar = PPatchPar;
    CppIsVec = [1 0 0 0 1 1 1];
    PModule = 'Valve';
    CppModule = 'Valve';
    P=loadPmoduledata(P,CA,PModule,CppModule,CppLoc,CppLoc1,PPatchPar,CppPar,CppIsVec);
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