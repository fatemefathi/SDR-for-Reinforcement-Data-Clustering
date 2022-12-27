function a = EpsilonGreedyActionSelection(centernumber, Q, Epsilon)

    nA = size(Q,2);
    
    if rand < Epsilon
        a = randi(nA);
    else
        QMax = max(Q(centernumber,:));
        A = find(Q(centernumber,:) == QMax);
        a = A(randi(numel(A)));
    end
    
end