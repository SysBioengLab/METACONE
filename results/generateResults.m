%% RESULTS: SECTION 4.1 METHODS 3.2 
% Metrics must be obtained for each model, and for each modality.
% For each model, we will run the algorithm five times (n = 5)
% The results will be LATER stored in the 'cases' variable.
% Each run will be stored and saved in RESULTADOS
% Each iteration will calculate the basis for each model CCi
% A log file will be saved to keep track of the iterations and executions.
% MODELS NEED TO BE ALREADY LOADED in the Models variable with their Names.

% ======================================================
% Change the path to the working directory.
cd /Users/alvin/Doctorado/TESIS/Article2/
% initCobraToolbox(false)
% changeCobraSolver('gurobi','ALL');
% ======================================================

% Parameters:
alpha = 0.1;
eTol = 1e-7;
% modality = 'full';
modality = 'fast';
% Names = {'ecolicore'};
% Models = {bmodel};

% Every model must have all consumption exchanges open open.
for i = 1:numel(Models)
    model = Models{i};
    ExIds = findExcRxns(model); 
    Exchanges = model.rxns(ExIds);
    if ismember(i, [3 4 5])
        Exchanges = Exchanges(contains(Exchanges, 'EX_'));
    end
    model = changeRxnBounds(model, Exchanges, -1e3, 'l');
    Models{i} = model;
end

% m = numel(Models);
% n = 5*m;
% nreps = 1:n;
% nombres = Names(1:m);
% %reps = split(num2str(1:5));
% reps = num2cell(1:5)';
% [Ns, Rs] = ndgrid(nombres, reps);
% cases = table(nreps',categorical(reshape(Ns',numel(Ns),1)), reshape(Rs',numel(Rs),1),...
%     'VariableNames',{'no.','model','rep'});
% cases.runtime = zeros(n,1);
% cases.noConv = zeros(n,1);
% cases.rank = zeros(n,1);
% cases.diversity = zeros(n,1);


% We run the algorithm 5 times permodel.
%parpool
for i = 1:5
    prenombre = ['RESULTS1_metacone_' modality]; 
    sufijo = '.mat';
    nombre = [prenombre '_' num2str(i) sufijo];

    try
        load(nombre);
    catch
        t = cell(numel(Models),2);
        RESULTADOS1 = struct();
        RESULTADOS1.vectors = t;
        RESULTADOS1.algoritmo = modality;
    end
    
    for j = 1:numel(Names)
        model = Models{j};
        [CC, Res] = metaCone(model,'Alpha',alpha,...
            'Nullity',true,'eTol',eTol,'Modality',modality);
        RESULTADOS1.vectors{j,1} = sparse(CC);
        RESULTADOS1.vectors{j,2} = Res;
        save(nombre, 'RESULTADOS1');
        Message1 = ['Finished ' num2str(i) '-th Iter. (' Names{j} '). Saved as ' nombre];
        fid = fopen(fullfile(cd, 'YourLogFile.txt'), 'a');
        fprintf(fid, '%s: %s\n', datestr(now, 0), Message1);
        fclose(fid);
    end
    
end; clear t model Message1 fid

%% Table 1

% Models characteristics 
nummos = numel(Names);
no_mets = zeros(nummos,1); % NUMBER OF METABOLITES
no_rxns = zeros(nummos,1); % NUMBER OF REACTIONS
for i = 1:(nummos)
    model = Models{i};
    no_mets(i) = size(model.S,1);
    no_rxns(i) = size(model.S,2);
end
table1 = table(Names, no_mets, no_rxns);
clear no_mets no_rxns model

%% Table 2

% PERFORMANCE OF metaCone ==================================================
% An alpha of 0.1 will be considered for THE ALGORITHM

% We take the saved results
n = 5*numel(Names); nreps = 1:n; reps = (1:5)';
[Ns, Rs] = ndgrid(Names, reps);
cases = array2table(zeros(n,3+8), 'VariableNames',...
    split('nreps model rep runtime noConv rank partMets div avgcs medcs iqrcs'));
cases(:,[1 3]) = table(nreps',categorical(reshape(Rs',numel(Rs),1)));
cases.model = categorical(reshape(Ns',numel(Ns),1));

for i = nreps

    prenombre = ['RESULTS1_metacone_' modality];
    sufijo = '.mat';
    nombre = [prenombre '_' num2str(cases.rep(i)) sufijo];
    load(nombre);
    name = cases.model(i);
    nid = name == Names;
    CC = full(RESULTADOS1.vectors{nid,1});
    R = RESULTADOS1.vectors{nid,2};
    
    %Participating Metabolites And diversity
    pMs = sum(any(CC'));
    Div = size(unique(sign(CC)','rows'),1);
    [~,cossims] = cosineSimilarity_vectors(CC);
    
    % Saving results
    cases(i,4:end) = {R.runtime R.noconv rank(CC,1e-9) pMs Div ...
        mean(cossims) median(cossims) iqr(cossims)};
end; clear RESULTADOS1 name nid CC R pMs Div cossims prenombre sufijo nombre Ns Rs i n nreps
% sound(gongo.y, gongo.Fs)

%Properties of the found basis
table2 = groupsummary(cases,"model",["mean" "std"],4:size(cases,2));
table2 = table2([5 3 1 4 2 6],[1 3 4 5 7 9 11 13 15 17]);
% colnames = split('model runtime noConv rank div cosSim nullity no_exch');
% table2 = array2table(zeros(nummos, size(colnames,1)));

RESULTADOS1.tabla1 = table1;
RESULTADOS1.tabla2 = table2;


%% SUPPLEMENTARY FILE S1. FIGURE S1. 
% Load the models first, with the main file.

% load RESULTS1_metacone_fast_1.mat
load RESULTS1_metacone_full_1.mat

figure(1)
for i = 1:size(RESULTADOS1.vectors,1)
    CC = RESULTADOS1.vectors{i,1};
    [cosSimMatrix, cossim] = cosineSimilarity_vectors(CC);
    subplot(2,3,i)
    histogram(cossim)
    title(Names{i})
    hold on
    xlabel("Cosine(\theta)")
    ylabel("Absolute Frequency")
end
hold off

%% SUPPLEMENTARY FILE S1. FIGURE S2.  
% We will analyse the modified e_coli core model.

bmodel = Models{2};
[cc, r] =  metaCone(bmodel,...
    'Nullity',true,...
    'eTol',1e-7,...
    'Alpha',0.1,...
    'Modality',modality, ...
    'keepAll',true);
clc; clf;

% BLOCK A ---------------------------------------------------------------
figure(2)
subplot(2,2,1)
% We first calculate the "cumulative rank" and fetch the epsilons.
allSols = r.Allsol;
nas = size(allSols, 2);
cumRank = zeros(nas,1);
nEps = zeros(nas,1); 
for i = 1:nas
    % SELECTED TOLERANCE 1e-9
    cumRank(i) = rank(full(allSols(:,1:i)),1e-9); % cumulative rank
    nEps(i) = norm(r.Epsilons(:,i)); %norms of Epsilons
end

plot(1:nas, nEps','o','MarkerEdgeColor',[.4 .4 .4],'MarkerSize',6)
%title("Rank and norm(Epsilons)")
xlim([0 size(r.Allsol,2)+2]); xlabel("Found Solutions")
ylim([-1 max(nEps)+1]); ylabel("Norm of epsilons (projections)")
hold on
yyaxis right
plot(1:nas, cumRank,'Color',"#4d4d9f",'LineWidth',3) %rgba(77,77,159,255)
ylabel("Cumulative rank of sols.")
ylim([0 rank(cc)+1])
ax = gca;
ax.FontSize = 15;
hold off

% BLOCK B ---------------------------------------------------------------
subplot(2,2,2)
% We calculate cosine similarity.
C = cc'*cc;
normas = sqrt(diag(C));
denom = normas*normas';
cosSim = C./denom;
%we add format to the figure.
vbc = [77 77 159]; % color?
Blulet = [linspace(1,vbc(1)/255,64)' ...
    linspace(1,vbc(2)/255,64)' ...
    linspace(1,vbc(3)/255,64)'];
heatmap(cosSim, 'Colormap', Blulet); grid off; 
%title("Cosine similarity between found conversions ")
ax = gca;
ax.FontSize = 15;

% BLOCK C ---------------------------------------------------------------
figure(2); subplot(2,2,3)
nus = nas;
cosSimsK = zeros(nus,1);
normas2 = zeros(nus,1);
% We calculate cosine similarity between a conversion and the 'next'
% conversion in the C_all matrix, accoring to indeces.
% We also compute the norm of the epsilon vector corresponding to the
% conversion.
for i = 1:(nus-1)
    v1 = allSols(:,i);
    v2 = allSols(:,i+1);
    cosSimsK(i) = dot(v1,v2)/(norm(v1)*norm(v2));
    normas2(i) = norm(r.Epsilons(:,i));
end
normas2(end) = norm(r.Epsilons(:,end));
%size(r.Allsol,2)-sum(cosSimsK>0.6)
[U,S,V] = svd(full(r.Allsol));
singulares = diag(S);
log10sings = log10(singulares);

%We will only plot diversity with cosine sim.
%Urank = zeros(nus,1);
Uniqueness = zeros(nus,1);
Allsol = full(r.Allsol);
% We calculate diversity using the unique() function.
for i = 1:nus
    actual = Allsol(:,1:i);
    sAllsol = sign(actual);
    %Urank(i) = rank(sAllsol);
    Uniqueness(i) = size(unique(sAllsol','rows'),1);
end
% plotting
plot(1:nus, cosSimsK','*')
%title("Orange: Diversity | Blue: Cosine")
xlim([0 size(r.Allsol,2)+5])
ylim([-.15 1.05])
ylabel("Cosine between v_i and v_{i+1}")
hold on
%plot(1:size(S,2), log10(singulares))
yyaxis right
plot(1:nus, Uniqueness,'LineWidth',3,'Color','#4d4d9f')
ylabel('Nº of unique sign columns')
legend({'Cos','Diversity'},'Location','best')
%ylim([0 100])
%plot(1:nus, normas2','o--')
ax = gca;
ax.FontSize = 15;
hold off

% BLOCK D ---------------------------------------------------------------
figure(2); subplot(2,2,4)
%hold on
plot(1:length(singulares), log10(singulares),'LineWidth',3)
xlabel("Index of singular value")
ylabel("Log_{10} of singular value")
xlim([0 23])
ylim([-30 6])
%yyaxis right
%plot(1:nus, Urank)
ax = gca;
ax.FontSize = 15;
%hold off
