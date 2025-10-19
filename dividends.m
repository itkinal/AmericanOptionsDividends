function qq = dividends(t,rq, exDiscrPropDates, propAmouns, ...
    exDiscrCashDates, cashAmouns)

global useDiscretePropDiv useDiscriteCashDiv

qq = rq(t);  
if useDiscretePropDiv && ~isempty(exDiscrPropDates) && ~isempty(propAmouns)
    iExDiv = 1;
    for s=2:length(t)                
        if iExDiv > length(exDiscrPropDates)                
            break
        elseif t(s-1) < exDiscrPropDates(iExDiv) && t(s) > exDiscrPropDates(iExDiv)
            qq(s) = qq(s) + propAmouns(iExDiv);  
            iExDiv = iExDiv + 1;
        end
    end
end

end
