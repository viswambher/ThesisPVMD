function panel = createDiagShadeMat1(ROWS,COLS,DIV,NSTEPS,WIDTH,angle,irr_u,irr_s,varargin)

% Sintax:
% panel = createDiagShadeMat1(ROWS,COLS,DIV,NSTEPS,angle,irr_u,irr_s)
% panel = createDiagShadeMat1(ROWS,COLS,DIV,NSTEPS,angle,irr_u,irr_s,plotting)
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
% panel = createDiagShadeMat1(ROWS,COLS,DIV,NSTEPS,WIDTH,angle,irr_u,irr_s,true);

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



