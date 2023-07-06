function [Mdl,impure_leaves,num_nodes,learn_time] = FnBuildDT(X,Y)
%FnBuildDT trains a DT using feature vectors provided by X and classifiers
%provided by Y

    tic
    Mdl = fitctree(X,Y,'CategoricalPredictors','all','MinParentSize',1);
    learn_time = toc;
    
    %Get number of nodes and impure leaves
    impure_leaves = nnz(Mdl.NodeRisk(~Mdl.IsBranchNode));
    num_nodes = Mdl.NumNodes;
end

