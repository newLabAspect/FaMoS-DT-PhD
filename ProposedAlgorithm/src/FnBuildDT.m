function [Mdl,impure_leaves,num_nodes,learn_time] = FnBuildDT(X,Y)
%FNBUILDDT Summary of this function goes here
%   Detailed explanation goes here
    tic
    Mdl = fitctree(X,Y,'CategoricalPredictors','all','MinParentSize',1);
    learn_time = toc;
    
    %Get number of nodes and impure leaves
    impure_leaves = nnz(Mdl.NodeRisk(~Mdl.IsBranchNode));
    num_nodes = Mdl.NumNodes;
end

