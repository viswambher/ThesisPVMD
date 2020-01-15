function panel = createDiagShadeMat(ROWS,COLS,DIV,NSTEPS,WIDTH,angle,irr_u,irr_s,varargin)

% Sintax:
% panel = createDiagShadeMat(ROWS,COLS,DIV,NSTEPS,angle,irr_u,irr_s)
% panel = createDiagShadeMat(ROWS,COLS,DIV,NSTEPS,angle,irr_u,irr_s,plotting)
%
% Description:
% Explain the purpose of the function
%
% %Example of use:
% 
% ROWS=8;
% COLS=12;
% DIV=10;
% NSTEPS = 10;
% angle = 0;%(deg)
% irr_s = 100;
% irr_u = 1000;
% WIDTH = 1;
% panel = createDiagShadeMat(ROWS,COLS,DIV,NSTEPS,WIDTH,angle,irr_u,irr_s,true);

if length(varargin)==1 && islogical(varargin{1})
    plotting = varargin{1};
else
    plotting = false;
end

if angle==90
    % m = tand(angle/2)
    b_init = 0;
    b_end = -COLS*DIV; %- m*COLS*DIV;
else
    m = tand(angle);
    b_init = 0;
    b_end=-ROWS*DIV-m*COLS*DIV;
end

if WIDTH<=0 && angle ~= 0 && angle~= 90 % There is definitely a better way to write this
    bs = linspace(b_init,b_end,NSTEPS+2);
    bs = bs(2:end-1);
    ks = ones(length(bs),1)*b_end;
% elseif WIDTH<=0 && angle ==0
%     bs = linspace(b_init,b_end,NSTEPS-1); %check
%     bs = bs(1:end);                %check
%     ks = ones(length(bs),1)*b_end;
% elseif  WIDTH<=0 && angle ==90
%     bs = linspace(b_init,b_end,NSTEPS-1); %check
%     bs = bs(1:end);                %check
%     ks = ones(length(bs),1)*b_end;

else
    cs = linspace(b_init,b_end,NSTEPS+1);
    cs = cs(1:end-1) -(cs(1)-cs(2))/2;
    bs = cs + WIDTH*DIV/2;
    ks = cs - WIDTH*DIV/2;
end

irr_d = irr_u-irr_s;

panel = zeros(ROWS,COLS,NSTEPS);

% b_init
% b_end
% NSTEPS
% bs
% cs
% ks

for ix=1:NSTEPS
    b=bs(ix);
    k=ks(ix);
    pexp = zeros(DIV*ROWS,DIV*COLS);
    for x=1:COLS*DIV
        for y=1:ROWS*DIV
            if angle>=90 && WIDTH<=0
                pexp(y,x)= irr_s + (-x<b)*irr_d;
            elseif angle < 90 && WIDTH <= 0
                pexp(y,x) = irr_s+(-y<m*x+b)*irr_d; %irradiance at each division
            elseif angle < 90 && WIDTH > 0 % WORKING
                pexp(y,x) = irr_s + (-y>m*x+k && -y<m*x+b)*(irr_d);
            elseif angle >= 90 && WIDTH > 0
                pexp(y,x) = irr_s + (-x>=k && -x<b)*(irr_d);
            else 
                pexp(y,x) = irr_s + rand(8,12).*irr_d
            end
        end
    end
% for ix=1:NSTEPS
%     b=bs(ix);
%     k=ks(ix);
%     pexp = zeros(DIV*ROWS,DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%             pexp(y,x) = irr_u-(-y>=m*x+k && -y<m*x+b)*(irr_u-irr_s); %irradiance at each division
%         end
%     end
    
    for row=1:size(panel,1)
        for col=1:size(panel,2)
            aux = zeros(size(pexp));
            aux(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %aux=1 for every division
            aux = aux.*pexp;
            panel(row,col,ix)= sum(aux(:))/DIV^2;
        end
    end
    if plotting
        figure
        [X,Y,Z] = meshgrid(0:size(panel,2),0:size(panel,1),0);
        surf(X,Y,Z,panel(:,:,ix));view(0,-90);
        caxis([0 1000])
        colorbar
        axis equal
        axis tight
    end    
% end
end
end



% % % % % % % % % % % % % % % ROWS=8;
% % % % % % % % % % % % % % % COLS=12;
% % % % % % % % % % % % % % % DIV=10;
% % % % % % % % % % % % % % % NSTEPS = 10;
% % % % % % % % % % % % % % % angle = 45;%(deg)
% % % % % % % % % % % % % % % irr_s = 100;
% % % % % % % % % % % % % % % irr_u = 1000;
% % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % irr_pattrn = zeros(ROWS,COLS,DIV*2);
% % % % % % % % % % % % % % % ix=1;
% % % % % % % % % % % % % % % for width =1:1:2
% % % % % % % % % % % % % % %     panel = createDiagShadeMat(ROWS,COLS,DIV,NSTEPS,width,angle,irr_u,irr_s);
% % % % % % % % % % % % % % %     irr_pattrn(:,:,ix:ix+DIV-1) = panel;
% % % % % % % % % % % % % % %     ix=ix+DIV
% % % % % % % % % % % % % % % end
% % % % % % % % % % % % % % % 




%%%%% First Version %%%%%

% function panel = createDiagShadeMat(ROWS,COLS,DIV,angle,irr_u,irr_s,varargin)
% 
% % Sintax:
% % panel = createDiagShadeMat(ROWS,COLS,DIV,angle,irr_u,irr_s)
% % panel = createDiagShadeMat(ROWS,COLS,DIV,angle,irr_u,irr_s,plotting)
% %
% % Description:
% % Explain the purpose of the function
% %
% % % Example of use:
% % ROWS=8;
% % COLS=12;
% % DIV=10;
% % angle = 45;%(deg)
% % irr_s = 100;
% % irr_u = 1000;
% % panel = createDiagShadeMat(ROWS,COLS,DIV,angle,irr_u,irr_s,true)
% 
% if length(varargin)==1 && islogical(varargin{1})
%     plotting = varargin{1};
% else
%     plotting = false;
% end
% 
% m = tand(angle);
% b_init = 0;
% b_end=-COLS*DIV-m*ROWS*DIV;
% 
% NSTEPS = 10;
% bs = linspace(b_init,b_end,NSTEPS+2);
% bs = bs(2:end-1);
% 
% irr_d = irr_u-irr_s;
% 
% panel = zeros(ROWS,COLS,NSTEPS);
% 
% for ix=1:NSTEPS
%     b=bs(ix);
%     
%     pexp = zeros(DIV*ROWS,DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%             pexp(y,x) = irr_s+(-y>m*x+b)*irr_d; %irradiance at each division
%         end
%     end
%     
%     for row=1:size(panel,1)
%         for col=1:size(panel,2)
%             aux = zeros(size(pexp));
%             aux(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %aux=1 for every division
%             aux = aux.*pexp;
%             panel(row,col,ix)= sum(aux(:))/DIV^2;
%         end
%     end
%     if plotting
%         figure
%         [X,Y,Z] = meshgrid(0:size(panel,2),0:size(panel,1),0);
%         surf(X,Y,Z,panel(:,:,ix));view(0,-90);
%         caxis([0 1000])
%         colorbar
%         axis equal
%         axis tight
%     end    
% end
% 
% %function payel = createDiagShadeMat1(ROWS,COLS,DIV,angle,irr_u,irr_s,varargin)
% 
% payel = zeros(ROWS,COLS,NSTEPS);
% 
% for ix=1:NSTEPS
%     b=bs(ix);
%     
%     peyp = zeros(DIV*ROWS,DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%            % pexp(y,x) = irr_s+(-y>m*x+b)*irr_d; %irradiance at each division
%             peyp(y,x) = irr_s+(-y<m*x+b)*irr_d;
%         end
%     end
%     
%     for row=1:size(panel,1)
%         for col=1:size(panel,2)
%             auy = zeros(size(peyp));
%             auy(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %auy=1 for every division
%             auy = auy.*peyp;
%             payel(row,col,ix)= sum(auy(:))/DIV^2;
%         end
%     end
% 
% if plotting
%          figure
%         [X,Y,Z] = meshgrid(0:size(payel,2),0:size(payel,1),0);
%         surf(X,Y,Z,payel(:,:,ix));view(0,-90);
%         caxis([0 1000])
%         colorbar
%         axis equal
%         axis tight
% end
% end
% 
% %% Shade moving in X-direction
% 
% % ROWS =8
% % COLS = 12
% CSTEPS = 12;
% pacel = zeros(ROWS,COLS,CSTEPS);
% 
% c_init = 0;
% c_end=120;
% 
% 
% cs = linspace(c_init,c_end,CSTEPS+1);
% cs = cs(1:end);
% 
% for ix=1:CSTEPS
%     c=cs(ix);
%     
%     pecp = zeros(DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%            % pexp(y,x) = irr_s+(-y>m*x+b)*irr_d; %irradiance at each division
%             pecp(y,x) = irr_s+(x>c)*irr_d;
%         end
%     end
%     
%     for row=1:size(pacel,1)
%         for col=1:size(pacel,2)
%             auc = zeros(size(pecp));
%             auc(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %auy=1 for every division
%             auc = auc.*pecp;
%             pacel(row,col,ix)= sum(auc(:))/DIV^2;
%         end
%     end
% 
% if plotting
%          figure
%         [X,Y,Z] = meshgrid(0:size(pacel,2),0:size(pacel,1),0);
%         surf(X,Y,Z,pacel(:,:,ix));view(0,-90);
%         caxis([0 1000])
%         colorbar
%         axis equal
%         axis tight
% end
% end
% 
% 
% 
% %% Shade moving in Y-direction
% % ROWS = 8;
% % COLS = 12;
% RSTEPS = 8;
% parel = zeros(ROWS,COLS,RSTEPS);
% 
% r_init = 0;
% r_end=80;
% 
% 
% rs = linspace(r_init,r_end,RSTEPS+1);
% %rs = rs(1:end);
% 
% for ix=1:RSTEPS
%     r=rs(ix);
%     
%     perp = zeros(DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%            % pexp(y,x) = irr_s+(-y>m*x+b)*irr_d; %irradiance at each division
%             perp(y,x) = irr_s+(y>r)*irr_d;
%         end
%     end
%     
%     for row=1:size(parel,1)
%         for col=1:size(parel,2)
%             aur = zeros(size(perp));
%             aur(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %auy=1 for every division
%             aur = aur.*perp;
%             parel(row,col,ix)= sum(aur(:))/DIV^2;
%         end
%     end
% 
% if plotting
%          figure
%         [X,Y,Z] = meshgrid(0:size(parel,2),0:size(parel,1),0);
%         surf(X,Y,Z,parel(:,:,ix));view(0,-90);
%         caxis([0 1000])
%         colorbar
%         axis equal
%         axis tight
% end
% end
% 
% %% Block Shade with size of one cell moving in Y-direction
% 
% % DIV = 10;
% % ROWS = 8;
% % COLS = 12;
% RRSTEPS = 8 ;
% parrel = zeros(ROWS,COLS,RRSTEPS);
% 
% rr_init = 0;
% rr_end=80;
% 
% %RRSTEPS = 8;
% rs = linspace(rr_init,rr_end,RRSTEPS+1);
% rs = rs(1:end);
% rm = linspace(rr_init+ DIV ,rr_end + DIV,RRSTEPS+1);
% rm = rm(1:end);
% 
% for ix=1:RSTEPS
%     r=rs(ix);
%     s = rm(ix);
%     perrp = zeros(DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%            % pexp(y,x) = irr_s+(-y>m*x+b)*irr_d; %irradiance at each division
%             perrp(y,x) = irr_s+(y>r && y<s)*irr_d;
%         end
%     end
%     
%     for row=1:size(parrel,1)
%         for col=1:size(parrel,2)
%             aurr = zeros(size(perrp));
%             aurr(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %auy=1 for every division
%             aurr = aurr.*perrp;
%             parrel(row,col,ix)= sum(aurr(:))/DIV^2;
%         end
%     end
% 
% if plotting
%          figure
%         [X,Y,Z] = meshgrid(0:size(parrel,2),0:size(parrel,1),0);
%         surf(X,Y,Z,parrel(:,:,ix));view(0,-90);
%         caxis([0 1000])
%         colorbar
%         axis equal
%         axis tight
% end
% end
% 
% %% Block Shade with size of one cell moving in X-direction
% 
% % DIV = 10;
% % ROWS = 8;
% % COLS = 12;
% CCSTEPS = 12 ;
% paccel = zeros(ROWS,COLS,CCSTEPS);
% 
% cc_init = 0;
% cc_end=120;
% 
% %CCSTEPS = 12;
% rs = linspace(cc_init,cc_end,CCSTEPS+1);
% rs = rs(1:end);
% cm = linspace(cc_init+ DIV ,cc_end + DIV,CCSTEPS+1);
% cm = cm(1:end);
% 
% for ix=1:CCSTEPS
%     r=rs(ix);
%     t = cm(ix);
%     peccp = zeros(DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%            % pexp(y,x) = irr_s+(-y>m*x+b)*irr_d; %irradiance at each division
%             peccp(y,x) = irr_s+(x>r && x<t)*irr_d;
%         end
%     end
%     
%     for row=1:size(paccel,1)
%         for col=1:size(paccel,2)
%             aucc = zeros(size(peccp));
%             aucc(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %auy=1 for every division
%             aucc = aucc.*peccp;
%             paccel(row,col,ix)= sum(aucc(:))/DIV^2;
%         end
%     end
% 
% if plotting
%          figure
%         [X,Y,Z] = meshgrid(0:size(paccel,2),0:size(paccel,1),0);
%         surf(X,Y,Z,paccel(:,:,ix));view(0,-90);
%         caxis([0 1000])
%         colorbar
%         axis equal
%         axis tight
% end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Block Shade with size of two cells moving in Y-direction
% 
% % DIV = 10;
% % ROWS = 8;
% % COLS = 12;
% RRSTEPS = 8 ;
% parrel = zeros(ROWS,COLS,RRSTEPS);
% 
% rr_init = 0;
% rr_end=80;
% 
% %RRSTEPS = 8;
% rs = linspace(rr_init,rr_end,RRSTEPS+1);
% rs = rs(1:end);
% rm = linspace(rr_init+ 2*DIV ,rr_end + 2*DIV,RRSTEPS+1);
% rm = rm(1:end);
% 
% for ix=1:RSTEPS
%     r=rs(ix);
%     s = rm(ix);
%     perrp = zeros(DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%            % pexp(y,x) = irr_s+(-y>m*x+b)*irr_d; %irradiance at each division
%             perrp(y,x) = irr_s+(y>r && y<s)*irr_d;
%         end
%     end
%     
%     for row=1:size(parrel,1)
%         for col=1:size(parrel,2)
%             aurr = zeros(size(perrp));
%             aurr(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %auy=1 for every division
%             aurr = aurr.*perrp;
%             parrel(row,col,ix)= sum(aurr(:))/DIV^2;
%         end
%     end
% 
% if plotting
%          figure
%         [X,Y,Z] = meshgrid(0:size(parrel,2),0:size(parrel,1),0);
%         surf(X,Y,Z,parrel(:,:,ix));view(0,-90);
%         caxis([0 1000])
%         colorbar
%         axis equal
%         axis tight
% end
% end
% 
% %% Block Shade with size of two cells moving in X-direction
% 
% % DIV = 10;
% % ROWS = 8;
% % COLS = 12;
% CCSTEPS = 12 ;
% paccel = zeros(ROWS,COLS,CCSTEPS);
% 
% cc_init = 0;
% cc_end=120;
% 
% %CCSTEPS = 12;
% rs = linspace(cc_init,cc_end,CCSTEPS+1);
% rs = rs(1:end);
% cm = linspace(cc_init+ 2*DIV ,cc_end + 2*DIV,CCSTEPS+1);
% cm = cm(1:end);
% 
% for ix=1:CCSTEPS
%     r=rs(ix);
%     t = cm(ix);
%     peccp = zeros(DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%            % pexp(y,x) = irr_s+(-y>m*x+b)*irr_d; %irradiance at each division
%             peccp(y,x) = irr_s+(x>r && x<t)*irr_d;
%         end
%     end
%     
%     for row=1:size(paccel,1)
%         for col=1:size(paccel,2)
%             aucc = zeros(size(peccp));
%             aucc(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %auy=1 for every division
%             aucc = aucc.*peccp;
%             paccel(row,col,ix)= sum(aucc(:))/DIV^2;
%         end
%     end
% 
% if plotting
%          figure
%         [X,Y,Z] = meshgrid(0:size(paccel,2),0:size(paccel,1),0);
%         surf(X,Y,Z,paccel(:,:,ix));view(0,-90);
%         caxis([0 1000])
%         colorbar
%         axis equal
%         axis tight
% end
% end
% 
% %% Block Shade with size of one cell moving in X-Y direction
% ROWS=8;
% COLS=12;
% DIV=10;
% angle = 45;%(deg)
% irr_s = 600;
% irr_u = 1000;
% 
% % if length(varargin)==1 && islogical(varargin{1})
% %     plotting = varargin{1};
% % else
% %     plotting = false;
% % end
% 
% 
% m = tand(angle);
% b_init = 0;
% b_end=-COLS*DIV-m*ROWS*DIV;
% 
% 
% NSTEPS = 10;
% bs = linspace(b_init,b_end,NSTEPS+2);
% bs = bs(2:end-1);
% ks = bs - (mean(bs)/2.5);
% 
% irr_d = irr_u-irr_s;
% 
% panel = zeros(ROWS,COLS,NSTEPS);
% 
% for ix=1:NSTEPS
%     b=bs(ix);
%     k=ks(ix);
%     pexp = zeros(DIV*ROWS,DIV*COLS);
%     for x=1:COLS*DIV
%         for y=1:ROWS*DIV
%             pexp(y,x) = irr_u-(-y>m*x+b && -y<m*x+k)*irr_s; %irradiance at each division
%         end
%     end
%     
%     for row=1:size(panel,1)
%         for col=1:size(panel,2)
%             aux = zeros(size(pexp));
%             aux(1+(row-1)*DIV:DIV+(row-1)*DIV,1+(col-1)*DIV:DIV+(col-1)*DIV)=1; %aux=1 for every division
%             aux = aux.*pexp;
%             panel(row,col,ix)= sum(aux(:))/DIV^2;
%         end
%     end
%    % if plotting
%         figure
%         [X,Y,Z] = meshgrid(0:size(panel,2),0:size(panel,1),0);
%         surf(X,Y,Z,panel(:,:,ix));view(0,-90);
%         caxis([0 1000])
%         colorbar
%         axis equal
%         axis tight
%   %  end    
% end
% 
% 
% 
% 
% 
% 
% 
