function f =  parseInferFunc(k)
    switch(k)
        case 'loopy' 
            f = @UGM_Infer_LBP;
        case 'chain'
            f = @UGM_Infer_Chain;
        case 'tree'
            f = @UGM_Infer_Tree;
        case 'mf'
            f = @UGM_Infer_MeanField;
        otherwise
              error('learn_parameters:parseInferFunc', ...
                  strcat('Unknown infer function :', k));
    end
end









