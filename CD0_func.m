function CD0 = CD0_func(S)
    
    % Different estimates of CD0
    % CD0 = 0.0447;
    % CD0 = 0.03951 * ((S).^(-0.1));
    CD0 = 0.03 * (S.^(-0.1));

end