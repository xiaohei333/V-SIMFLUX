function result = checkPhase(pdata)
    result = mod(real(pdata+pi), 2*pi)-pi;
end