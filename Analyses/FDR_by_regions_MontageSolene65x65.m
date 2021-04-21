clear;clc
load('C:\Users\laura\AnalysesNIRS\OtherProjects\Solene\Facteurs_nana_sansmat_beauxsujets_2021_03_31_groupe_pval.mat')
NC=65;
channels=1:NC;
alpha = 0.05; %%%you can modify this%%%

letter= {'A' 'B' 'M' 'E' 'C' 'D'  'O' 'F' };
label={'F_R' 'FT_R' 'C_R' 'T_R'  'F_L' 'FT_L' 'C_L' 'T_L'  };
hemis={'right H' 'right H' 'right H' 'right H' 'right L' 'right L' 'right L' 'right L'};
roi={'frontal' 'frontotemporal' 'centra' 'temporal' 'frontal' 'frontotemporal' 'centra' 'temporal' };
antpos={'ant' 'ant' 'pos'  'pos' 'ant' 'ant' 'pos' 'pos'};

%make new matrice without the short-distance (P detector)
for c=1:NC
    if strcmp('P',ZoneList{c}(1))
        SD= c;
    end
end
channels(SD)=[];

%organize channel number according to ROI
for r=1:length(letter)
    roiChan(r).label=label{r};
    roiChan(r).hemis=hemis(r);
    roiChan(r).roi=roi(r);
    roiChan(r).antpos=antpos(r);
    roiChan(r).nchan=[];
    for c=channels
        if strcmp(letter{r},ZoneList{c}(1))
            roiChan(r).nchan=[roiChan(r).nchan c];
        end
    end
end

mat(1).label='all'; mat(1).listchan=ZoneList([roiChan(1:8).nchan]);
mat(1).pval=matcorr([roiChan(1:8).nchan],[roiChan(1:8).nchan]);

mat(2).label='left H'; mat(2).listchan=ZoneList([roiChan(5:8).nchan]);
mat(2).pval=matcorr([roiChan(5:8).nchan],[roiChan(5:8).nchan]);
mat(3).label='right H'; mat(3).listchan=ZoneList([roiChan(1:4).nchan]);
mat(3).pval=matcorr([roiChan(1:4).nchan],[roiChan(1:4).nchan]);

mat(4).label='anterior L_R'; mat(4).listchan=ZoneList([roiChan([1 2 5 6]).nchan]);
mat(4).pval=matcorr([roiChan([1 2 5 6]).nchan],[roiChan([1 2 5 6]).nchan]);
mat(5).label='posterior L_R';mat(5).listchan=ZoneList([roiChan([3 4 7 8]).nchan]);
mat(5).pval=matcorr([roiChan([3 4 7 8]).nchan],[roiChan([3 4 7 8]).nchan]);

mat(6).label='frontal L_R';mat(6).listchan=ZoneList([roiChan([1 5]).nchan]);
mat(6).pval=matcorr([roiChan([1 5]).nchan],[roiChan([1 5]).nchan]);
mat(7).label='frontotemporal L_R';mat(7).listchan=ZoneList([roiChan([2 6]).nchan]);
mat(7).pval=matcorr([roiChan([2 6]).nchan],[roiChan([2 6]).nchan]);
mat(8).label='central L_R';mat(8).listchan=ZoneList([roiChan([3 7]).nchan]);
mat(8).pval=matcorr([roiChan([3 7]).nchan],[roiChan([3 7]).nchan]);
mat(9).label='temporal L_R';mat(9).listchan=ZoneList([roiChan([4 8]).nchan]);
mat(9).pval=matcorr([roiChan([4 8]).nchan],[roiChan([4 8]).nchan]);

mat(10).label='anterior L';mat(10).listchan=ZoneList([roiChan([5 6]).nchan]);
mat(10).pval=matcorr([roiChan([5 6]).nchan],[roiChan([5 6]).nchan]);
mat(11).label='anterior R';mat(11).listchan=ZoneList([roiChan([1 2]).nchan]);
mat(11).pval=matcorr([roiChan([1 2]).nchan],[roiChan([1 2]).nchan]);
mat(12).label='posterior L';mat(12).listchan=ZoneList([roiChan([7 8]).nchan]);
mat(12).pval=matcorr([roiChan([7 8]).nchan],[roiChan([7 8]).nchan]);
mat(13).label='posterior R';mat(13).listchan=ZoneList([roiChan([3 4]).nchan]);
mat(13).pval=matcorr([roiChan([3 4]).nchan],[roiChan([3 4]).nchan]);

for m=1:length(mat)
    
    original=mat(m).pval;
    temp=[];
    for c=1:size(original,1)
        for cc=(c+1):size(original,2)
            temp=[temp original(c,cc)];
        end
    end
    
    %%%%fdr%%%%%%%%%%%%%
    
  
    % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
    [~, ~, ~, q1] = fdr_bh(temp, alpha, 'pdep', 'yes' );
    n_sig = sum(q1 <= alpha);
    mat(m).corrP{1,1}=n_sig;
    mat(m).corrP{1,2}=sprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini,1995)\n',n_sig, alpha);
    
    %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1... doesn't make sense
    [q2, ~] = bonf_holm(temp, alpha);
    n_sig = sum(q2 <= alpha);
    mat(m).corrP{2,1}=n_sig;
    mat(m).corrP{2,2}=sprintf('%d tests are significant p<=%.2f using Bonferroni-Holm correction\n',n_sig, alpha);
    
    %Multicmp function to apply Holm Step Down Procedure%%
    [q3,tmpalpha] = multicmp (temp','down',alpha);
    n_sig = sum(q3 <= alpha);
    mat(m).corrP{3,1}=n_sig;
    mat(m).corrP{3,2}=sprintf('%d tests are significant p<=%.2f using Holm Step Down Procedure\n',n_sig, alpha);
    
    %Multicmp function to apply Hochberg's step up procedure%%
    [q4,tmpalpha] = multicmp (temp','up',alpha);
    n_sig = sum(q4 <= alpha);
    mat(m).corrP{4,1}=n_sig;
    mat(m).corrP{4,2}=sprintf('%d tests are significant p<=%.2f using Hochberg Step Up Procedure\n',n_sig, alpha);
    
    %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
    [q5,tmpalpha] = multicmp (temp','fdr',0.05);
    n_sig = sum(q5 <= alpha);
    mat(m).corrP{5,1}=n_sig;
    mat(m).corrP{5,2}=sprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini, 1995)\n',n_sig, alpha);
    
      % % FDR according to Storey 2002%%%
    %DANS LE SCRIPT MAFDR... CA SEMBLE ÊTRE À LA LIGNE 212 QUE LES Q VALUES
    %SONT CALCULÉES... JE NE COMPRENDS PAS LE CALCUL, MAIS ÇA SEMBLE FAIRE
    %QQCH D'ÉTRANGE. 
%     [q6,test] = mafdr(temp);
%     n_sig = sum(q6 <= alpha);
%     mat(m).corrP{6,1}=n_sig;
%     mat(m).corrP{6,2}=sprintf('%d tests are significant p<=%.2f using FDR correction (Storey, 2002)\n',n_sig, alpha);
%     
    
    %reshape the matrices to be saved
    mat(m).corrP{1,3}=nan(size(original));
    mat(m).corrP{2,3}=nan(size(original));
    mat(m).corrP{3,3}=nan(size(original));
    mat(m).corrP{4,3}=nan(size(original));
    mat(m).corrP{5,3}=nan(size(original));

    xx=1;
    for c=1:size(original,1)
        for cc=(c+1):size(original,2)
            mat(m).corrP{1,3}(c,cc)=q1(xx);
            mat(m).corrP{2,3}(c,cc)=q2(xx);
            mat(m).corrP{3,3}(c,cc)=q3(xx);
            mat(m).corrP{4,3}(c,cc)=q4(xx);
            mat(m).corrP{5,3}(c,cc)=q5(xx);

            
           mat(m).corrP{1,3}(cc,c)=q1(xx);
            mat(m).corrP{2,3}(cc,c)=q2(xx);
            mat(m).corrP{3,3}(cc,c)=q3(xx);
            mat(m).corrP{4,3}(cc,c)=q4(xx);
            mat(m).corrP{5,3}(cc,c)=q5(xx);

            xx=xx+1;
        end
    end
    

[~,pos]=max([mat(m).corrP{1:5,1}]);
    mat(m).results=mat(m).corrP{pos,2};
    mat(m).adj_pval=mat(m).corrP{pos,3};
    clear q*
end

save('resultats_par_regions_2021-04-21.mat','roiChan','mat','letter','label','matcorr','ZoneList')

