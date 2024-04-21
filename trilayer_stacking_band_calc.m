clear all;

abc_stack = true; % false for ABA, true for ABC
BZ_mesh = 1; % 0 for high-sym line cut, 1 for BZ 2D mesh sampling

lattice_a=1.42*sqrt(3);

layer_d=[0,0,3.35];

% lattice_a*
a1=lattice_a*[sqrt(3)/2,-1/2]';
a2=lattice_a*[sqrt(3)/2,+1/2]';

A = [a1 a2];
G = 2*pi*inv(A)';

b2 = G(:,1)+G(:,2);
b1 = G(:,1);

K_pt = 1/3*(2*G(:,1)+G(:,2));

b_shift = (a1+a2)/3;

pos_list = zeros(3,6);

% B orbital positions
pos_list(1:2,2) = b_shift;
pos_list(1:2,4) = b_shift;
pos_list(1:2,6) = b_shift;

% AB stacking (first two layers)
pos_list(:,3) = pos_list(:,3) + [b_shift;layer_d(3)];
pos_list(:,4) = pos_list(:,4) + [b_shift;layer_d(3)];

% top layer
pos_list(:,5) = pos_list(:,5) + [0;0;2*layer_d(3)];
pos_list(:,6) = pos_list(:,6) + [0;0;2*layer_d(3)];

% ABC stacking
if (abc_stack)
    pos_list(:,5) = pos_list(:,5) + 2*[b_shift;0];
    pos_list(:,6) = pos_list(:,6) + 2*[b_shift;0];   
end

sc_grid = 5;

% graphene TBH from DFT with strain corrections
t_list = [-2.892, 0.243, -0.266]; % turn off 2nd term for p-h symmetric dirac cone
r_list = lattice_a/sqrt(3)*[1,sqrt(3),2];


if (BZ_mesh == 1)
    k_samp = 40;

    % Do a reduced sampling around the K point because full BZ wastes a lot of
    % time on high energy states.

    dmesh = 0.10*linspace(-1,1,3*k_samp+1);
    dmesh = dmesh(1:end-1);
    [mesh_x,mesh_y] = meshgrid(dmesh,dmesh);

    k_list(1,:) = K_pt(1)+b1(1)*mesh_x(:) + b2(1)*mesh_y(:);
    k_list(2,:) = K_pt(2)+b1(2)*mesh_x(:) + b2(2)*mesh_y(:);
else
    G_pt = 0.0*b2;
    M_pt = 0.5*b2;
    lc_points = [G_pt, K_pt, M_pt, G_pt];
    k_samp = 50;
    dmesh = linspace(0,1,k_samp+1);
    dmesh = dmesh(1:end-1);
    k_list = zeros(2,3*k_samp);
    for seg = 0:2
        start_idx = k_samp*seg;
        pt1 = lc_points(:,seg+1);
        pt2 = lc_points(:,seg+2);
        k_list(:,start_idx+[1:k_samp]) = pt1*(1-dmesh) + pt2*dmesh;
    end
    k_ax(1) = 0;
    for idx = 2:size(k_list,2)
        k_ax(idx) = k_ax(idx-1) + norm(k_list(:,idx) - k_list(:,idx-1));
    end
    
    
    
end

for k_idx = 1:size(k_list,2)

    kh = k_list(:,k_idx);
    H = zeros(6,6);
            
    % interlayer terms
    for i = -sc_grid:sc_grid
        for j = -sc_grid:sc_grid
            R = i*a1 + j*a2;
            for l = 1:2
                for o_from = 1:2
                    idx_f = (l-1)*2+o_from;
                    pos_f = pos_list(:,idx_f);
                    for o_to = 1:2
                        idx_t = (l)*2+o_to;
                        pos_t = pos_list(:,idx_t);

                        dr = [R;0] + pos_t - pos_f;
                        t_h = dft_interlayer_coupling(dr', 0, 0, lattice_a);
                        %t_h = 0; % turns off inter coupling
                        phase_h = exp(1j*dot(kh,dr(1:2)));
                        H(idx_t,idx_f) = H(idx_t,idx_f) + t_h*phase_h;

                    end

                end
            end
        end
    end
    H = H + H'; % get hermitian conj. of interlayer terms
            
            % monolayer terms
    for i = -sc_grid:sc_grid
        for j = -sc_grid:sc_grid
            R = i*a1 + j*a2;
            for l = 1:3
                for o_from = 1:2
                    idx_f = (l-1)*2+o_from;
                    pos_f = pos_list(1:2,idx_f);
                    for o_to = 1:2
                        idx_t = (l-1)*2+o_to;
                        pos_t = pos_list(1:2,idx_t);

                        dr = R + pos_t - pos_f;
                        rl = norm(dr);

                        t_h = 0;
                        for idx = 1:3
                            if abs(rl - r_list(idx)) < 1e-4
                               t_h = t_list(idx);
                            end
                        end

                        phase_h = exp(1j*dot(kh,dr(1:2)));
                        H(idx_t,idx_f) = H(idx_t,idx_f) + t_h*phase_h;

                    end

                end
            end
                       

        end
    end

    bands(k_idx,:) = sort(real(eigs(H,6)));
    
end

%%
if (BZ_mesh == 1)
    E_f = -0.729; %  to move Fermi energy back to 0 eV
    [dos_sweep, ~, E_list, ~] = interp_kp_dos([1], {bands'-E_f}, {k_list'});
end
%%

clf

figure;
subplot(2,1,1)
if (BZ_mesh == 0)
    plot(k_ax,bands-E_f,'k')
    axis([0,max(k_ax),-11,11])
    set(gca,'XTick',[0, k_ax(1+k_samp),k_ax(2+2*k_samp)])
    xticklabels({'G','K','M'})
else
    plot(bands,'k')
end
%axis([80 82 -2 2])
if (abc_stack)
    title('ABC trilayer graphene bands (near K)')    
else
    title('ABA trilayer graphene bands (near K)')
end
ylabel('Energy (eV)')
xlabel('k index')

if (BZ_mesh == 1)
    % spin and valley already taken care of in interp_kp_dos.m!
    dos = dos_sweep{1};
    

    subplot(2,1,1)
    plot(bands-E_f,'-k','Linewidth',0.1)
    axis([1 inf -0.7 0.7])
    if (abc_stack)
        title('ABC trilayer graphene BZ sampling')    
    else
        title('ABA trilayer graphene BZ sampling')
    end
    ylabel('Energy (eV)')
    
    subplot(2,1,2)
    plot(E_list,dos,'k')
    %axis([-0.5 0.5 0 1.5])
    if (abc_stack)
        title('ABC trilayer graphene DOS')    
    else
        title('ABA trilayer graphene DOS')
    end
    xlabel('Energy (eV)')
    ylabel('DOS (states per (eV nm$^2$))')
    set(gca,'Xtick',[-0.5:0.1:0.5])
end

%% Save variables

if (abc_stack)
    save('ABC_highsymbands.mat','bands', 'k_samp', 'E_f')
    save('ABC_dos.mat', 'E_list', 'dos')
else
    save('ABA_highsymbands.mat','bands', 'k_samp', 'E_f')
    save('ABA_dos.mat', 'E_list', 'dos')
end

%% Load test

load('ABC_highsymbands.mat')
figure;
plot(k_ax,bands,'k')
axis([0,max(k_ax),-11,11])
set(gca,'XTick',[0, k_ax(1+k_samp),k_ax(2+2*k_samp)])
xticklabels({'G','K','M'})
ylabel('Energy (eV)')
xlabel('k index')

%% Plot dos

load('ABC_dos.mat')
figure;
plot(E_list,dos,'k')
title('ABC trilayer graphene DOS')
xlabel('Energy (eV)')
ylabel('DOS (states per (eV nm$^2$))', 'interpreter','latex')
set(gca,'Xtick',[-0.5:0.1:0.5])
xlim([-0.5, 0.5])