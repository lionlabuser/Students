classdef GlobalRemover
    % Global component remove program is made by Dr. Zhang, Xian,  Brain function lab, Yale school of medicine, Psychiatry department
    % Please site my paper: "Separation of the global and local components in functional near-infrared spectroscopy signals using principal component spatial filtering"
    % and "Signal processing of functional NIRS data acquired during overt speaking" which indicating global component in dexoyHb signal
    % xzhang63@gmail.com
    properties
        distance_matrix
        sigma
        kernel
        vraw
        polhemus
        xyz
        nCh
        badCh
    end % end of propertieso
    methods
        function x = GlobalRemover(polhemus,sigma,badCh)
            x.polhemus = polhemus;
            x.badCh=badCh;
            x.sigma=sigma;
            ind=(1+polhemus.nTR):size(polhemus.NFRI_result.OtherC,1);
            badCh(badCh>length(ind))=[];
            ind(badCh)=[];
            x.xyz=polhemus.NFRI_result.OtherC(ind,:);
            x.nCh=size(x.xyz,1);
            [TH,PHI,R] = cart2sph(x.xyz(:,1),x.xyz(:,2),x.xyz(:,3));
            thphi(:,1)=TH;
            thphi(:,2)=PHI;
            x.distance_matrix =squareform(pdist(thphi,@arclen));
            kernel=exp(- x.distance_matrix.^2 / x.sigma);
            for i=1:x.nCh
                x.kernel(:,i)=kernel(:,i)/sum(kernel(:,i));
            end
        end
        function v= getGlobal(x,vraw)
            for i=1:x.nCh
                v(i)=vraw(1:x.nCh)'*x.kernel(:,i);
            end
        end
        function v=remove(x,vraw)
            v=vraw'-x.getGlobal(vraw);
            v=v';
        end
    end
end
