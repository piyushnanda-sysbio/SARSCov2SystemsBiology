expr=readtable('D:/COVID19/Data/FPKM.csv')

cd D:/COVID19/Model/GSMs
load 'Recon3DOut.mat'
model
ids=readtable('GeneExtName.csv')
model.genes=table2array(ids(:,2))
model.genes

model=changeRxnBounds(model,model.rxns(findExcRxns(model)),0,'l')
unique(model.lb(findExcRxns(model)))

%ChangeCOBRASolver

changeCobraSolver('ibm_cplex')

%Change bounds as per the media composition (HAM)
cd ..
media=readtable('HAMcompositionInHMR.txt')
exc=table2array(media(:,2))
exc=exc(~cellfun('isempty',exc))
model=changeRxnBounds(model,exc,-0.00133,'l')
model=changeObjective(model,'biomass_reaction',1)
sol=optimizeCbModel(model)

names=table2array(expr(:,1))
met_genes=names(ismember(names,model.genes))

expressionData_UI.gene=string(met_genes)
expressionData_I.gene=string(met_genes)

geneExpr_UI=table2array(expr(:,2:4))
geneExpr_UI=geneExpr_UI(ismember(table2array(expr(:,1)),met_genes))

geneExpr_I=table2array(expr(:,5:7))
geneExpr_I=geneExpr_I(ismember(table2array(expr(:,1)),met_genes))

val_UI=mean(geneExpr_UI,2)
val_I=mean(geneExpr_I,2)
val_UI=(val_UI-mean(val_UI,1))/std(val_UI)
val_I=(val_I-mean(val_I,1))/std(val_I)
expressionData_UI.value=val_UI
expressionData_I.value=val_I

hist(expressionData_I.value,50)
xlim([min(expressionData_I.value), max(expressionData_I.value)])
hist(expressionData_UI.value,50)
xlim([min(expressionData_UI.value), max(expressionData_UI.value)])

[expressionRxns_UI parsedGPR_UI] = mapExpressionToReactions(model, expressionData_UI)
[expressionRxns_I parsedGPR_I] = mapExpressionToReactions(model, expressionData_I)

hist(expressionRxns_I,50)
xlim([min(expressionRxns_I), max(expressionRxns_I)])
hist(expressionRxns_UI,50)
xlim([min(expressionRxns_UI), max(expressionRxns_UI)])

p = 0:0.25:1;
y1 = quantile(expressionRxns_I,p);
y2 = quantile(expressionRxns_UI,p);
z_i = [p;y1]
z_ui = [p;y2]

changeCobraSolver('ibm_cplex','milp')

tissueUI = iMAT(model, expressionRxns_UI, -0.2, 0, 1E-8, {}, 'MILPLog', 10800, 1)
tissueI = iMAT(model, expressionRxns_I, -0.2, 0, 1E-8, {}, 'MILPLog', 10800, 1)

rxn_for=printRxnFormula(model,'biomass_reaction')
tissueUI=addReaction(tissueUI,'biomass_reaction','reactionFormula',char(rxn_for))
tissueI=addReaction(tissueI,'biomass_reaction','reactionFormula',char(rxn_for))

tissueUI=changeObjective(tissueUI,'biomass_reaction',1)
tissueI=changeObjective(tissueI,'biomass_reaction',1)

BiomassComponents = model.mets(find(model.S(:, strmatch('biomass_reaction', model.rxns))))
[tissueUI_NEW, rxnNames] = addDemandReaction(tissueUI, BiomassComponents);
for i = 1 : length(BiomassComponents)
        tissueUI_NEW = changeObjective(tissueUI_NEW, rxnNames(i),1)
        FBAsolution = optimizeCbModel(tissueUI_NEW, 'max')
        BiomassComponentsValueUI(i,1) = FBAsolution.f
end

BiomassComponents = model.mets(find(model.S(:, strmatch('biomass_reaction', model.rxns))))
[tissueI_NEW, rxnNames] = addDemandReaction(tissueI, BiomassComponents);
for i = 1 : length(BiomassComponents)
        tissueI_NEW = changeObjective(tissueI_NEW, rxnNames(i),1)
        FBAsolution = optimizeCbModel(tissueI_NEW, 'max')
        BiomassComponentsValueI(i,1) = FBAsolution.f
end
table(BiomassComponents, BiomassComponentsValueI,BiomassComponentsValueUI)


save 'tissueUI.mat'
save 'tissueI.mat'


