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
    roiChan(r).updateChanMAT=((r-1)*8+1):(8*r);
end

%% SECTION TARGETING THE SAME CHANNELS IN THE VERTICAL AND HORIZZONTAL AXES
mat(1).label='all'; mat(1).listchanV=ZoneList([roiChan(1:8).nchan]);
mat(1).pval=matcorr([roiChan(1:8).nchan],[roiChan(1:8).nchan]);

mat(2).label='left H'; mat(2).listchanV=ZoneList([roiChan(5:8).nchan]);
mat(2).pval=matcorr([roiChan(5:8).nchan],[roiChan(5:8).nchan]);
mat(3).label='right H'; mat(3).listchanV=ZoneList([roiChan(1:4).nchan]);
mat(3).pval=matcorr([roiChan(1:4).nchan],[roiChan(1:4).nchan]);

mat(4).label='anterior L_R'; mat(4).listchanV=ZoneList([roiChan([1 2 5 6]).nchan]);
mat(4).pval=matcorr([roiChan([1 2 5 6]).nchan],[roiChan([1 2 5 6]).nchan]);
mat(5).label='posterior L_R';mat(5).listchanV=ZoneList([roiChan([3 4 7 8]).nchan]);
mat(5).pval=matcorr([roiChan([3 4 7 8]).nchan],[roiChan([3 4 7 8]).nchan]);

mat(6).label='frontal L_R';mat(6).listchanV=ZoneList([roiChan([1 5]).nchan]);
mat(6).pval=matcorr([roiChan([1 5]).nchan],[roiChan([1 5]).nchan]);
mat(7).label='frontotemporal L_R';mat(7).listchanV=ZoneList([roiChan([2 6]).nchan]);
mat(7).pval=matcorr([roiChan([2 6]).nchan],[roiChan([2 6]).nchan]);
mat(8).label='central L_R';mat(8).listchanV=ZoneList([roiChan([3 7]).nchan]);
mat(8).pval=matcorr([roiChan([3 7]).nchan],[roiChan([3 7]).nchan]);
mat(9).label='temporal L_R';mat(9).listchanV=ZoneList([roiChan([4 8]).nchan]);
mat(9).pval=matcorr([roiChan([4 8]).nchan],[roiChan([4 8]).nchan]);

mat(10).label='anterior L';mat(10).listchanV=ZoneList([roiChan([5 6]).nchan]);
mat(10).pval=matcorr([roiChan([5 6]).nchan],[roiChan([5 6]).nchan]);
mat(11).label='anterior R';mat(11).listchanV=ZoneList([roiChan([1 2]).nchan]);
mat(11).pval=matcorr([roiChan([1 2]).nchan],[roiChan([1 2]).nchan]);
mat(12).label='posterior L';mat(12).listchanV=ZoneList([roiChan([7 8]).nchan]);
mat(12).pval=matcorr([roiChan([7 8]).nchan],[roiChan([7 8]).nchan]);
mat(13).label='posterior R';mat(13).listchanV=ZoneList([roiChan([3 4]).nchan]);
mat(13).pval=matcorr([roiChan([3 4]).nchan],[roiChan([3 4]).nchan]);

mat(14).label='intraH only frontal R';mat(14).listchanV=ZoneList([roiChan([1]).nchan]);
mat(14).pval=matcorr([roiChan([1]).nchan],[roiChan([1]).nchan]);
mat(15).label='intraH only frontotemporal R';mat(15).listchanV=ZoneList([roiChan([2]).nchan]);
mat(15).pval=matcorr([roiChan([2]).nchan],[roiChan([2]).nchan]);
mat(16).label='intraH only central R';mat(16).listchanV=ZoneList([roiChan([3]).nchan]);
mat(16).pval=matcorr([roiChan([3]).nchan],[roiChan([3]).nchan]);
mat(17).label='intraH only temporal R';mat(17).listchanV=ZoneList([roiChan([4]).nchan]);
mat(17).pval=matcorr([roiChan([4]).nchan],[roiChan([4]).nchan]);

mat(18).label='intraH only frontal L';mat(18).listchanV=ZoneList([roiChan([5]).nchan]);
mat(18).pval=matcorr([roiChan([5]).nchan],[roiChan([5]).nchan]);
mat(19).label='intraH only frontotemporal L';mat(19).listchanV=ZoneList([roiChan([6]).nchan]);
mat(19).pval=matcorr([roiChan([6]).nchan],[roiChan([6]).nchan]);
mat(20).label='intraH only central L';mat(20).listchanV=ZoneList([roiChan([7]).nchan]);
mat(20).pval=matcorr([roiChan([7]).nchan],[roiChan([7]).nchan]);
mat(21).label='intraH only temporal L';mat(21).listchanV=ZoneList([roiChan([8]).nchan]);
mat(21).pval=matcorr([roiChan([8]).nchan],[roiChan([8]).nchan]);
for m=1:length(mat)
    
    %copy vertical channel list to horizontal channel list
    mat(m).listchanH=mat(m).listchanV;
    
    
    %reshape temporary the matrice for 1D instead of 2D
    original=mat(m).pval;
    temp=[];
    for c=1:size(original,1)
        for cc=(c+1):size(original,2)
            temp=[temp original(c,cc)];
        end
    end
    
    %calculate number of sig results WItHOUT correction
    n_sig = sum(temp <= alpha);
     mat(m).res0=sprintf('%d tests are significant p<=%.2f WITHOUT correction\n',n_sig, alpha);
    
    
    %%%%fdr%%%%%%%%%%%%%
    
    
    % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
    [~, ~, ~, q1] = fdr_bh(temp, alpha, 'pdep', 'yes' );
    n_sig = sum(q1 <= alpha);
    
    mat(m).res1=sprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini,1995)\n',n_sig, alpha);
    
    %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1... doesn't make sense
    [q2, ~] = bonf_holm(temp, alpha);
    n_sig = sum(q2 <= alpha);
    
    mat(m).res2=sprintf('%d tests are significant p<=%.2f using Bonferroni-Holm correction\n',n_sig, alpha);
    
    %Multicmp function to apply Holm Step Down Procedure%%
    [q3,tmpalpha] = multicmp (temp','down',alpha);
    n_sig = sum(q3 <= alpha);
    
    mat(m).res3=sprintf('%d tests are significant p<=%.2f using Holm Step Down Procedure\n',n_sig, alpha);
    
    %Multicmp function to apply Hochberg's step up procedure%%
    [q4,tmpalpha] = multicmp (temp','up',alpha);
    n_sig = sum(q4 <= alpha);
    
    mat(m).res4=sprintf('%d tests are significant p<=%.2f using Hochberg Step Up Procedure\n',n_sig, alpha);
    
    %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
    [q5,tmpalpha] = multicmp (temp','fdr',0.05);
    n_sig = sum(q5 <= alpha);
    
    mat(m).res5=sprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini, 1995)\n',n_sig, alpha);
    
    % % FDR according to Storey 2002%%%
    %DANS LE SCRIPT MAFDR... CA SEMBLE ÊTRE À LA LIGNE 212 QUE LES Q VALUES
    %SONT CALCULÉES... JE NE COMPRENDS PAS LE CALCUL, MAIS ÇA SEMBLE FAIRE
    %QQCH D'ÉTRANGE.
    %     [q6,test] = mafdr(temp);
    %     n_sig = sum(q6 <= alpha);
    %     mat(m).corrP{6,1}=n_sig;
    %     mat(m).corrP{6,2}=sprintf('%d tests are significant p<=%.2f using FDR correction (Storey, 2002)\n',n_sig, alpha);
    %
    
    
    
    %reshape the matrice from the temporary 1D into 2D
    %start by empty matrices with nan
    mat(m).corrP1=nan(size(original));
    mat(m).corrP2=nan(size(original));
    mat(m).corrP3=nan(size(original));
    mat(m).corrP4=nan(size(original));
    mat(m).corrP5=nan(size(original));
    
    xx=1;
    for c=1:size(original,1)
        for cc=(c+1):size(original,2)
            mat(m).corrP1(c,cc)=q1(xx);
            mat(m).corrP2(c,cc)=q2(xx);
            mat(m).corrP3(c,cc)=q3(xx);
            mat(m).corrP4(c,cc)=q4(xx);
            mat(m).corrP5(c,cc)=q5(xx);
            
            
            mat(m).corrP1(cc,c)=q1(xx);
            mat(m).corrP2(cc,c)=q2(xx);
            mat(m).corrP3(cc,c)=q3(xx);
            mat(m).corrP4(cc,c)=q4(xx);
            mat(m).corrP5(cc,c)=q5(xx);
            
            xx=xx+1;
        end
    end
    
    clear q*
end



%% SECTION TARGETING DIFFERENT CHANNELS IN THE VERTICAL AXIS AND HORIZZONTAL AXIS
z=22;
mat(z).label='interH all';
targetvertical=[roiChan([1:4]).nchan];
targethorizon=[roiChan([5:8]).nchan];
mat(z).listchanV=ZoneList(targetvertical); mat(z).listchanH=ZoneList(targethorizon);
mat(z).pval=matcorr(targetvertical,targethorizon);

z=23;
mat(z).label='interH anterior';
targetvertical=[roiChan([1 2]).nchan];
targethorizon=[roiChan([5 6]).nchan];
mat(z).listchanV=ZoneList(targetvertical); mat(z).listchanH=ZoneList(targethorizon);
mat(z).pval=matcorr(targetvertical,targethorizon);

z=24;
mat(z).label='interH posterior';
targetvertical=[roiChan([3 4]).nchan];
targethorizon=[roiChan([7 8]).nchan];
mat(z).listchanV=ZoneList(targetvertical); mat(z).listchanH=ZoneList(targethorizon);
mat(z).pval=matcorr(targetvertical,targethorizon);

z=25;
mat(z).label='interH only frontal L_H';
targetvertical=[roiChan([1]).nchan];
targethorizon=[roiChan([5]).nchan];
mat(z).listchanV=ZoneList(targetvertical); mat(z).listchanH=ZoneList(targethorizon);
mat(z).pval=matcorr(targetvertical,targethorizon);

z=26;
mat(z).label='interH only frontotemporal L_H';
targetvertical=[roiChan([2]).nchan];
targethorizon=[roiChan([6]).nchan];
mat(z).listchanV=ZoneList(targetvertical); mat(z).listchanH=ZoneList(targethorizon);
mat(z).pval=matcorr(targetvertical,targethorizon);

z=27;
mat(z).label='interH only central L_H';
targetvertical=[roiChan([3]).nchan];
targethorizon=[roiChan([7]).nchan];
mat(z).listchanV=ZoneList(targetvertical); mat(z).listchanH=ZoneList(targethorizon);
mat(z).pval=matcorr(targetvertical,targethorizon);

z=28;
mat(z).label='interH only temporal L_H';
targetvertical=[roiChan([4]).nchan];
targethorizon=[roiChan([8]).nchan];
mat(z).listchanV=ZoneList(targetvertical); mat(z).listchanH=ZoneList(targethorizon);
mat(z).pval=matcorr(targetvertical,targethorizon);

for m=22:length(mat)
    
    
    %reshape temporary the matrice for 1D instead of 2D - but we keep ALL
    %VALUES
    original=mat(m).pval;
    temp=reshape(original,[1 numel(original)]);
    
     %calculate number of sig results WItHOUT correction
    n_sig = sum(temp <= alpha);
     mat(m).res0=sprintf('%d tests are significant p<=%.2f WITHOUT correction\n',n_sig, alpha);
    
    %%%%fdr%%%%%%%%%%%%%
    
    % FDR according to Benjamini & Hochberg(1995) de Groppe%%%%
    [~, ~, ~, q1] = fdr_bh(temp, alpha, 'pdep', 'yes' );
    n_sig = sum(q1 <= alpha,'all');
    
    mat(m).res1=sprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini,1995)\n',n_sig, alpha);
    
    %Bonferroni-Holm correction by Groppe%%% Corrected_p are greater than 1... doesn't make sense
    [q2, ~] = bonf_holm(temp, alpha);
    n_sig = sum(q2 <= alpha,'all');
    
    mat(m).res2=sprintf('%d tests are significant p<=%.2f using Bonferroni-Holm correction\n',n_sig, alpha);
    
    %Multicmp function to apply Holm Step Down Procedure%%
    [q3,tmpalpha] = multicmp (temp','down',alpha);
    n_sig = sum(q3 <= alpha,'all');
    
    mat(m).res3=sprintf('%d tests are significant p<=%.2f using Holm Step Down Procedure\n',n_sig, alpha);
    
    %Multicmp function to apply Hochberg's step up procedure%%
    [q4,tmpalpha] = multicmp (temp','up',alpha);
    n_sig = sum(q4 <= alpha,'all');
    mat(m).corrP4=n_sig;
    mat(m).res4=sprintf('%d tests are significant p<=%.2f using Hochberg Step Up Procedure\n',n_sig, alpha);
    
    %Multicmp function to apply FDR correction according to Benjamini & Hochberg%%
    [q5,tmpalpha] = multicmp (temp','fdr',0.05);
    n_sig = sum(q5 <= alpha,'all');
    
    mat(m).res5=sprintf('%d tests are significant p<=%.2f using FDR correction (Benjamini, 1995)\n',n_sig, alpha);
    
    % % FDR according to Storey 2002%%%
    %DANS LE SCRIPT MAFDR... CA SEMBLE ÊTRE À LA LIGNE 212 QUE LES Q VALUES
    %SONT CALCULÉES... JE NE COMPRENDS PAS LE CALCUL, MAIS ÇA SEMBLE FAIRE
    %QQCH D'ÉTRANGE.
    %     [q6,test] = mafdr(temp);
    %     n_sig = sum(q6 <= alpha);
    %     mat(m).corrP{6,1}=n_sig;
    %     mat(m).corrP{6,2}=sprintf('%d tests are significant p<=%.2f using FDR correction (Storey, 2002)\n',n_sig, alpha);
    %
    
    
    
    %reshape the matrice from the temporary 1D into 2D
    
    mat(m).corrP1=reshape(q1,size(original));
    mat(m).corrP2=reshape(q2,size(original));
    mat(m).corrP3=reshape(q3,size(original));
    mat(m).corrP4=reshape(q4,size(original));
    mat(m).corrP5=reshape(q5,size(original));
    
    
    clear q*
end


%% SAVE RESULTS
save(['resultats_par_regions_' datestr(now,'yyyy-mm-dd') '.mat'],'roiChan','mat','letter','label','matcorr','ZoneList')

