function [DB, J, J_DB, X] = start_from_pre_results(oldname, n_bins)

    old_mat = load(oldname);
    DB = old_mat.DB;
    J = old_mat.J;
    J_DB = old_mat.J_DB;
    X = DB(2).X;
    X(1:3,n_bins+1)=0;
