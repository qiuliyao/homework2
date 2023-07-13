
slice=zeros(N-1);
nt_max=200;%max no. of iteration 
er=zeros(nt_max,1);

     b=1.6; %optimun value of the SOR method, 
    onemb=1-b;
    
    nt=0; %iteration number
    while nt<nt_max
        nt=nt+1;
        
        for jj=1:N
            
            if(jj==1)
                jjm=N; %jj-1
            else
                jjm=jj-1;
            end
            
            if(jj==N)
                jjp=1;
            else
                jjp=jj+1;
            end
            
            
            for kk=1:N
                
                if(kk==1)
                    kkm=N;
                else
                    kkm=kk-1;
                end
                if(kk==N)
                    kkp=1;
                else
                    kkp=kk+1;
                end
                
                % interior points
                psi(kk,jj)=0.25*b*(psi(kkm,jj)+psi_guess(kkp,jj)+psi(kk,jjm)+psi_guess(kk,jjp)-dx*dx*zeta0(kk,jj))    +   .....
                    onemb*psi_guess(kk,jj);
                
            end
        end



        %diagnose the error with respect the zeta_true
 if mod( nt,50)==0
        slice = (psi(ix+1,jy) +psi(ix-1,jy) +psi(ix,jy+1) +psi(ix,jy-1)- 4*psi(ix,jy))*r_dx2;
        slice = slice-zeta0(2:end-1,2:end-1);
        
        slice=log10(abs( slice));
        er(nt)=max(slice(:));
         if er(nt)<-12
%              disp([num2str(nt),' SOR iters'])
            break
        end
end



        psi_guess=psi;
    end
    
% % % %     for diagnostics    % % % %
%     figure;
%     plot(er,'b.')
%     title('log10(error in zeta0)')
    
%     figure;
%     imagesc(slice);colorbar
%     title('log10(diff error)')










 












