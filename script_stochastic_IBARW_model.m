%--------------------------------------------------------------------------
% Inflationary theory of branching morphogenesis in the mouse salivary gland
% Bordeu I, Chatzeli L, and Simons BD (2022).
%--------------------------------------------------------------------------
% Stochastic simulation of the three-dimensional (3d) inflationary branch-
% ing-arresting random walk (IBARW), for uniform and isotropic tissue
% expansion.
%
% See supplementary Note for details.
%
% INPUTS: see MODEL PARAMETERS section below, and Methods in paper. The
% default parameters provided correspond to the parameters used to simulate
% the E14.5 to E18.5 development of the SG using one of three E14.5
% rudimentary trees as initial conditions
%
% OUTPUTS:
% edge_list: Array of size (number_of_links*3), containing the source node
%            id (column 1), target node id (column 2), and distance from
%            source to target (column 3).
% node_positions : Array of size number_of_nodes*4. Column 1: node id.
%                  Columns 2-4: coordinates (x,y,z) for each node, respectively.
% edge_list_tree: Reduced edge_list, where each node corresponds to either
%                 a branching point or a termination point or the root node
% node_positions_tree : node positions for the reduced branching tree.
%
% For questions: ib443 (at) cam.ac.uk
%--------------------------------------------------------------------------
clear all;close all;clc;

%% Working directory
% define a working directory to save output files (only used if save_files == 1)
% in the form 'D:\Documents\BARW_model\sim_output\'
work_path = 'E:\Documents\Ignacio\Salivary gland\scripts-Salivary-gland_git\simulations\For publication\';
save_files = 0; % set to 1 if want to save output files in work_path folder
plot_over_time = 1; % set to 1 if want to plot simulation (set n_reps to 1)

%% SIMULATION PARAMETERS -------------------------------------------------
boundary_type = 'open'; 
n_reps = 1; % number of realisations

% MODEL PARAMETERS
sigma = 0.05; % R_branch: branching ratio
expansion_rate = [0.006]; % R_exp: rate of expansion (h^-1)
annihil_radius = 45; % R_a: annihilation radius (h^-1)
t_max = [84]; % maximum simuation time (h)
% seed types: 
%    'tree E14.5': initializes the system with E14.5 rudimentary trees
%    'single seed': initializes the system with a single active tip
seed_type = 'tree E14.5'; 
% Other parameters
v0 = 1; % um/h
noise_amplitude = pi/7; % sqrt(2*D_r), with D_r the rotational.diff.const.
sensing_angle = pi/3; % angle at which tips snese ducts
dt = 1;
%--------------------------------------------------------------------------
dL = dt*v0;

% Distribution of branch angles (alpha):
pretruncMean = 90; pretruncSD = 20;
untruncated = makedist('Normal',pretruncMean,pretruncSD);
truncated = truncate(untruncated,50,150); % used in simulations

% Distribution of roll angles (phi): mean and SD (sqrt(2)b) for the Laplace distirbution 
mean_lap = 90;
std_lap = sqrt(2)*21;

for n_rep = 1:n_reps
    main_lengths = [];
    [init_active_part,tot_part_number,active_part_number,active_part_pos,active_part_angle,edge_list,node_positions,edge_list_clean,node_positions_clean] = rudimentary_tree(seed_type,n_reps,n_rep);

    init_inactive_part = 0;
    inactive_part_number = [];
    inactive_part_pos = [];
    inactive_part_angle = [];

    init_arrested_part = 0;
    arrested_part_number = [];
    arrested_part_pos = [];
    inactive_part_angle = [];
    %%
    t = 0;

    n_active_part = init_active_part;
    n_inactive_part = init_inactive_part;

    new_part_pos = [];

    while t < t_max && n_active_part>0
        % update time
        t = t + dt;
        % expand system
        if expansion_rate > 0

            [active_part_pos,inactive_part_pos,edge_list,node_positions,edge_list_clean,node_positions_clean,tot_part_number] = expand_domain(active_part_pos,inactive_part_pos,edge_list,node_positions,edge_list_clean,node_positions_clean,expansion_rate,tot_part_number,dL,dt);

            % reactivation of inactive tips due to tissue expansion
            [active_part_pos,active_part_angle,active_part_number,n_active_part, inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part] = check_reactivated_particles(active_part_pos,active_part_angle,active_part_number,n_active_part, inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,annihil_radius,node_positions,sensing_angle);
        end

        % random ordering of tip particles
        part_order = randperm(n_active_part);
        % draw random number fto perform
        % actions
        ran = rand(n_active_part,1);
        active_part_number = active_part_number(part_order);
        active_part_pos = active_part_pos(part_order,:);
        active_part_angle = active_part_angle(part_order,:);

        for part_num = n_active_part:-1:1
            % choose a tip cell to elongate or branch

            if isempty(edge_list) || ran(part_num) <= 1/(1+sigma) % tip elongation
                % find location and angle of tip to branch
                part_pos = active_part_pos(part_num,2:end);
                part_angle = active_part_angle(part_num,2:end);
                % compute new location of the tip
                [new_part_pos, new_part_angle] = update_particle_position(part_pos,part_angle,dL,dt,noise_amplitude);
                % check distance of particle to ducts
                h = check_distance_to_trace(part_pos,new_part_pos,annihil_radius,node_positions,sensing_angle);
                % update location in list
                [edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,tot_part_number] = add_particle(h,edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,tot_part_number,part_num,new_part_pos,new_part_angle);
            else % branch
                % find location and angle of tip to branch
                part_pos = active_part_pos(part_num,2:end);
                part_angle = active_part_angle(part_num,2:end);

                % branching angle from truncated normal distribution
                ang = random(truncated)*pi/180;
                branch_angles = [ang/2 -ang/2];
                % roll angle from laplace distribution
                rotate_angle = laprnd(1, 1, mean_lap, std_lap)*pi/180;

                % compute position of two new tips
                [new_part_pos1,new_part_angle1,new_part_pos2,new_part_angle2] = branch_tip(part_pos,part_angle,rotate_angle,branch_angles,dL);
                % check if any of the two new tips is close to a duct
                h1 = check_distance_to_trace(part_pos,new_part_pos1,annihil_radius,node_positions,sensing_angle);
                h2 = check_distance_to_trace(part_pos,new_part_pos2,annihil_radius,node_positions,sensing_angle);
                % add new tips to list
                [edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,tot_part_number] = add_two_particles(h1,h2,edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,tot_part_number,part_num,new_part_pos1,new_part_angle1,new_part_pos2,new_part_angle2,dL);

            end
        end

        if plot_over_time == 1
            subplot(1,2,1)
            % plot network
            plot_graph_3d(edge_list,node_positions);
            hold on
            % add actie and inactive tips
            if n_inactive_part > 0
                scatter3(inactive_part_pos(:,2),inactive_part_pos(:,3),inactive_part_pos(:,4),10,'b','filled')
            end
            if n_active_part > 0
                scatter3(active_part_pos(:,2),active_part_pos(:,3),active_part_pos(:,4),10,'r','filled')
            end
            hold off

            axis equal; box on
            set(gcf,'color','w'); view(2)
            title(['t = ',num2str(t),' h'])
            xlabel('\mum');ylabel('\mum');zlabel('\mum')
            
            subplot(1,2,2)
            source_node = 1;
            [edge_list_large,node_positions_large] = clean_edge_list(edge_list_clean,node_positions_clean);
            G = graph(edge_list_large(:,1),edge_list_large(:,2));
            h = plot(G,'-','NodeLabel',{},'Marker','none','EdgeColor','k');
            layout(h,'layered','sources',source_node);
            ylabel('level i');
            box on
            pause(0.001)
        end
    end
end

% -------------------------------------------------------------------------
% FUNCTIONS ---------------------------------------------------------------
% -------------------------------------------------------------------------

function [new_pos,new_angle] = update_particle_position(part_pos,part_angle,dL,dt,noise_amplitude)
part_angle(1) = part_angle(1) + 2*(rand()-0.5)*sqrt(dt)*noise_amplitude;
part_angle(2) = part_angle(2) + (acos(1 - 2*rand())/2-1/2)*sqrt(dt)*noise_amplitude;
new_pos = part_pos + dL*[sin(part_angle(1))*cos(part_angle(2)) sin(part_angle(1))*sin(part_angle(2)) cos(part_angle(1))];
new_angle = part_angle;
end

function [new_pos1,new_angle1,new_pos2,new_angle2] = branch_tip(part_pos,part_angle,rotate_angle,branch_angles,dL)

% normalised reference axis:
ref_ax = [sin(part_angle(1))*cos(part_angle(2)) sin(part_angle(1))*sin(part_angle(2)) cos(part_angle(1))];
% random rotation angle
%     rotate_angle = pi*rand();
% branch angles
theta1 = part_angle(1)+branch_angles(1);
theta2 = part_angle(1)+branch_angles(2);
% generate new branches
pos = dL*[sin(theta1)*cos(part_angle(2)) sin(theta1)*sin(part_angle(2)) cos(theta1)];
v1 = rotate_3D(pos', 'any', rotate_angle, ref_ax')';
pos = dL*[sin(theta2)*cos(part_angle(2)) sin(theta2)*sin(part_angle(2)) cos(theta2)];
v2 = rotate_3D(pos', 'any', rotate_angle, ref_ax')';

new_angle1 = [atan2(sqrt(v1(1)^2 + v1(2)^2),v1(3)), atan2(v1(2),v1(1))];
new_angle2 = [atan2(sqrt(v2(1)^2 + v2(2)^2),v2(3)), atan2(v2(2),v2(1))];

new_pos1 = v1 + part_pos;
new_pos2 = v2 + part_pos;
end

function [edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,tot_part_number] = add_particle(h,edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,tot_part_number,part_num,new_part_pos,new_part_angle)
% Add new particle postions to the list and keep as inactive if close to a
% duct (h = 1)
if h
    inactive_part_pos(end+1,:) = active_part_pos(part_num,:);
    inactive_part_angle(end+1,:) = active_part_angle(part_num,:);
    inactive_part_number(end+1) = active_part_number(part_num);
    n_inactive_part = n_inactive_part + 1;

    active_part_pos(part_num,:) = [];
    active_part_angle(part_num,:) = [];
    active_part_number(part_num) = [];
    n_active_part = n_active_part - 1;
else
    tot_part_number = tot_part_number+1;
    edge_list(end+1,:) = [active_part_pos(part_num,1) tot_part_number 1];

    [r,c] = find(edge_list_clean == active_part_pos(part_num,1));
    if numel(r) == 0
        edge_list_clean(end+1,:) = [active_part_pos(part_num,1) tot_part_number 1];
        node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos];
    else
        edge_list_clean(r,c) = tot_part_number;
        edge_list_clean(r,3) = edge_list_clean(r,3) + 1;

        node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos];
    end

    active_part_pos(part_num,:) = [tot_part_number,new_part_pos];
    active_part_angle(part_num,:) = [tot_part_number,new_part_angle];
    node_positions(end+1,:) = [tot_part_number new_part_pos];
    active_part_number(part_num) = tot_part_number;
end
end

function [edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,tot_part_number] = add_two_particles(h1,h2,edge_list,node_positions,edge_list_clean,node_positions_clean,active_part_pos,active_part_angle,active_part_number,n_active_part,inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,tot_part_number,part_num,new_part_pos1,new_part_angle1,new_part_pos2,new_part_angle2,dL)
% Here we add both new tip particles to the list
edge_list(end+1,:) = [active_part_pos(part_num,1) tot_part_number+1 dL];
edge_list(end+1,:) = [active_part_pos(part_num,1) tot_part_number+2 dL];
edge_list_clean(end+1,:) = [active_part_pos(part_num,1) tot_part_number+1 dL];
edge_list_clean(end+1,:) = [active_part_pos(part_num,1) tot_part_number+2 dL];

if h1==0 && h2==0
    % If both particles are far from duct, then keep them as active
    tot_part_number = tot_part_number+1;
    active_part_pos(part_num,:) = [tot_part_number,new_part_pos1];
    active_part_angle(part_num,:) = [tot_part_number,new_part_angle1];
    node_positions(tot_part_number,:) = [tot_part_number new_part_pos1];
    node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos1];
    active_part_number(part_num) = tot_part_number;

    tot_part_number = tot_part_number+1;
    active_part_pos(end+1,:) = [tot_part_number,new_part_pos2];
    active_part_angle(end+1,:) = [tot_part_number,new_part_angle2];
    node_positions(tot_part_number,:) = [tot_part_number new_part_pos2];
    node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos2];
    active_part_number(end+1) = tot_part_number;

    n_active_part = n_active_part + 1;

elseif h1==0 && h2==1
    % If one particle is arrested, then keep one as active
    tot_part_number = tot_part_number+1;
    active_part_pos(part_num,:) = [tot_part_number,new_part_pos1];
    active_part_angle(part_num,:) = [tot_part_number,new_part_angle1];
    node_positions(tot_part_number,:) = [tot_part_number new_part_pos1];
    node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos1];
    active_part_number(part_num) = tot_part_number;

    tot_part_number = tot_part_number+1;
    inactive_part_pos(end+1,:) = [tot_part_number,new_part_pos2];
    inactive_part_angle(end+1,:) = [tot_part_number,new_part_angle2];
    node_positions(tot_part_number,:) = [tot_part_number new_part_pos2];
    node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos2];
    inactive_part_number(end+1) = tot_part_number;

    n_inactive_part = n_inactive_part + 1;

elseif h1==1 && h2==0
    % If one particle is arrested, then keep one as active
    tot_part_number = tot_part_number+1;
    inactive_part_pos(end+1,:) = [tot_part_number,new_part_pos1];
    inactive_part_angle(end+1,:) = [tot_part_number,new_part_angle1];
    node_positions(tot_part_number,:) = [tot_part_number new_part_pos1];
    node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos1];
    inactive_part_number(end+1) = tot_part_number;

    n_inactive_part = n_inactive_part + 1;

    tot_part_number = tot_part_number+1;
    active_part_pos(part_num,:) = [tot_part_number,new_part_pos2];
    active_part_angle(part_num,:) = [tot_part_number,new_part_angle2];
    node_positions(tot_part_number,:) = [tot_part_number new_part_pos2];
    node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos2];
    active_part_number(part_num) = tot_part_number;


elseif h1==1 && h2==1
    % If both particles are arrested, keep both as inactive
    tot_part_number = tot_part_number+1;
    inactive_part_pos(end+1,:) = [tot_part_number,new_part_pos1];
    inactive_part_angle(end+1,:) = [tot_part_number,new_part_angle1];
    node_positions(tot_part_number,:) = [tot_part_number new_part_pos1];
    node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos1];
    inactive_part_number(end+1) = tot_part_number;

    n_inactive_part = n_inactive_part + 1;

    tot_part_number = tot_part_number+1;
    inactive_part_pos(end+1,:) = [tot_part_number,new_part_pos2];
    inactive_part_angle(end+1,:) = [tot_part_number,new_part_angle2];
    node_positions(tot_part_number,:) = [tot_part_number new_part_pos2];
    node_positions_clean(tot_part_number,:) = [tot_part_number new_part_pos2];
    inactive_part_number(end+1) = tot_part_number;

    active_part_pos(part_num,:) = [];
    active_part_angle(part_num,:) = [];
    active_part_number(part_num) = [];
    n_active_part = n_active_part - 1;

end
end

function [edge_list_large,node_positions_large] = clean_edge_list(edge_list_large,node_positions_large)
% this is to put the nodes in edge_list_clean and node_positions_clean in
% sequencial order.
un_nodes = unique(edge_list_large(:,1:2));
complete_list = 1:max(un_nodes);
rm_nodes = setdiff(complete_list,un_nodes); % these are the nodes that must be deleted
for i = 1:numel(un_nodes)
    node_positions_large(un_nodes(i),1) = i;
    ind = (edge_list_large(:,1) == un_nodes(i));
    edge_list_large(ind,1) = i;
    ind = (edge_list_large(:,2) == un_nodes(i));
    edge_list_large(ind,2) = i;
end
node_positions_large(rm_nodes,:) = [];

end % end clean_edge_list

function h = check_distance_to_trace(part_pos,new_part_pos,annihil_radius,node_positions,sensing_angle)
% Here we check the distance and angle from a tip to the ducts, using its
% propagation direction as vetor of reference.
all_pos = node_positions(:,2:end);
inds = (abs(all_pos(:,1)-new_part_pos(1)) < annihil_radius & abs(all_pos(:,2)-new_part_pos(2)) < annihil_radius & abs(all_pos(:,3)-new_part_pos(3)) < annihil_radius);

trace_pos = all_pos(inds,:);
dists = pt_to_pts_dist(new_part_pos,trace_pos);

v = trace_pos - new_part_pos;
w = new_part_pos - part_pos;
w = repmat(w,size(v,1),1);
angles = atan2(vecnorm(cross(v,w,2),2,2),dot(v,w,2))';

h = any(dists < annihil_radius & abs(angles) < sensing_angle, 'all');
end

function [init_active_part,tot_part_number,active_part_number,active_part_pos,active_part_angle,edge_list,node_positions,edge_list_clean,node_positions_clean] = rudimentary_tree(option,n_reps,n_rep)
% Here we construct the initial condition either with a single stalk
% ('single seed') or loading an E14.5 rudimentary tree ('tree E14.5').
switch option
    case 'single seed'
        init_pos = zeros(1,3);
        init_angle = zeros(1,2); % this angle indicates the direction of linear growth
        init_angle(1) = pi/2;
        init_active_part = 1;

        % active particle positions:
        active_part_pos = [1,init_pos];
        active_part_angle = [1,init_angle];
        % total number of particles
        tot_part_number = 1;
        % number of ative tips:
        active_part_number = 1;

        % Network info:
        edge_list = [];
        node_positions = [1, init_pos];

        edge_list_clean = [];
        node_positions_clean = [1, init_pos];

    case 'tree E14.5'
        if ~ispc; slsh = '/'; else; slsh = '\'; end

        if n_rep <= n_reps/3 % 1/3 of realisations are initialised with each one of three templates.
            main_folder_path = ['E14.5_trees',slsh','Sample 1',slsh];
            dr = [0.5679, 0.5679, 1]; % resolution of each experimentla sample (um/pixel)
        elseif n_rep <= 2*n_reps/3
            main_folder_path = ['E14.5_trees',slsh','Sample 2',slsh];
            dr = [0.5678, 0.5678, 3]; % (um/pixel)
        else
            main_folder_path = ['E14.5_trees',slsh','Sample 3',slsh];
            dr = [0.5679, 0.5679, 2]; % (um/pixel)
        end

        % load template
        load([main_folder_path,'edge_list_large.dat']);
        load([main_folder_path,'node_positions_large.dat']);

        % lets centre the template
        node_positions_large(:,2:4) = (node_positions_large(:,2:4)-node_positions_large(1,2:4)).*dr; % bring to origin and normalize

        n1 = edge_list_large(:,1);
        n2 = edge_list_large(:,2);
        edge_list_large(:,3) = sqrt(sum(dr.^2.*(node_positions_large(n1,2:4)-node_positions_large(n2,2:4)).^2,2));

        % lets put the endnodes in the last positions in the list
        [~,~,D] = node_level(edge_list_large,node_positions_large,1);
        ends = find(D == 1); ends(ends == 1) = []; % node 1 is the root
        M = max(node_positions_large(:,1));
        anodes = (M-length(ends)+1):M;
        anodes = setdiff(anodes,ends);
        node_pos_temp = node_positions_large;
        edge_list_temp = edge_list_large;
        for i = 1:length(anodes)
            edge_list_large(edge_list_temp == ends(i)) = anodes(i);
            edge_list_large(edge_list_temp == anodes(i)) = ends(i);
            node_positions_large(ends(i),2:4) = node_pos_temp(anodes(i),2:4);
            node_positions_large(anodes(i),2:4) = node_pos_temp(ends(i),2:4);
        end

        [~,~,D,G] = node_level(edge_list_large,node_positions_large,1);

        % extract active endnode info (angle, position, number...)
        active_part_angle = [];
        active_part_number = find(D == 1); active_part_number(active_part_number == 1) = [];
        init_active_part = length(active_part_number);
        active_part_pos = node_positions_large(active_part_number,:);
        % extract angles
        for i = 1:size(active_part_pos,1)
            apnum = active_part_pos(i,1);
            app = active_part_pos(i,2:4);

            [r,c] = find(edge_list_large(:,1:2) == apnum);
            if c == 1; c0 = 2; else c0 = 1; end

            apnum0 = edge_list_large(r,c0);
            app0 = node_positions_large(apnum0,2:4);

            u = app-app0;
            angle1 = atan2(u(2),u(1));
            angle2 = atan2(norm(u(1:2)),u(3));

            active_part_angle(end+1,:) = [apnum, angle2, angle1];
        end

        [edge_list,node_positions] = add_intermediate_nodes(edge_list_large,node_positions_large,1);

        tot_part_number = length(node_positions(:,1));

        edge_list_clean = edge_list_large;
        node_positions_clean = node_positions_large;
end
end % end rudimentary_tree

function [edge_list,node_positions] = add_intermediate_nodes(edge_list,node_positions,dr)
% Add intermediate nodes when distance between any two nodes is larger than
% 1 unit.
r = 0;
while ~isempty(r)
    r = find(edge_list(:,3) > dr);
    tot_part_number = max(node_positions(:,1));
    if ~isempty(r)
        n1 = edge_list(r,1);
        n2 = edge_list(r,2);
        d = edge_list(r,3);

        new_parts_num = [(tot_part_number + 1):(tot_part_number+length(r))]';
        tot_part_number = tot_part_number + length(r);
        edge_list(end+1:end+length(r),:) = [n1 new_parts_num d./2.0];
        edge_list(end+1:end+length(r),:) = [new_parts_num n2 d./2.0];
        pos_new = (node_positions(n1,2:4) + node_positions(n2,2:4))./2.0;
        node_positions(new_parts_num,:) = [new_parts_num pos_new];
        % if tot_part_number ~= length(

        edge_list(r,:) = [];
    end
end 
end % end add_intermediate_nodes

function [G,h] = plot_graph_3d(edge_list,node_positions)
% Plot 3d network
m = max(node_positions(:,1));
edge_list(any(isnan(edge_list), 2), :) = m + 1;
G = graph(edge_list(:,1),edge_list(:,2));
h = plot(G,'XData',node_positions(:,2),'YData',node_positions(:,3),'ZData',node_positions(:,4),'NodeLabel',{},'Marker','none','EdgeColor','k');
end % end plot_graph_3d

function [paths,dists,D,G] = node_level(edge_list,node_positions,source_node)
% Here we extract network information:
% paths: list of node connecting each node with the source_node
% dists: distance (in level) from each node to the source_node
% D: node degree list
% G: graph object
G = graph(edge_list(:,1),edge_list(:,2));
paths = {};
dists = [];
for i = 1:size(node_positions,1)
    nod = node_positions(i,1);
    if nod~=i
        disp('Warning: incompatible node_positions file')
        break;
    end
    if isfinite(nod)
        [P,d] = shortestpath(G,source_node,nod);
        paths{i} = P;
        dists(i) = d;
    else
        paths{i} = nan;
        dists(i) = nan;
    end
end
D = degree(G);
end % end node_level

function [active_part_pos,inactive_part_pos,edge_list,node_positions,edge_list_clean,node_positions_clean,tot_part_number] = expand_domain(active_part_pos,inactive_part_pos,edge_list,node_positions,edge_list_clean,node_positions_clean,expansion_rate,tot_part_number,dL,dt)
% here we expand the system by rescaling the location of all particles by (1+dt*expansion_rate)
if ~isempty(edge_list)
    node_positions(:,2:4) = node_positions(:,2:4).*(1+dt*expansion_rate);
    node_positions_clean(:,2:4) = node_positions_clean(:,2:4).*(1+dt*expansion_rate);
    active_part_pos(:,2:4) = active_part_pos(:,2:4).*(1+dt*expansion_rate);
    if ~isempty(inactive_part_pos)
        inactive_part_pos(:,2:4) = inactive_part_pos(:,2:4).*(1+dt*expansion_rate);
    end
    
    % add intermediate nodes if distance between nodes is larger than dL
    % after the expansion.
    n1 = edge_list(:,1);
    n2 = edge_list(:,2);
    edge_list(:,3) = sqrt(sum(dL^2.*(node_positions(n1,2:4)-node_positions(n2,2:4)).^2,2));
    n1 = edge_list_clean(:,1);
    n2 = edge_list_clean(:,2);
    edge_list_clean(:,3) = sqrt(sum(dL^2.*(node_positions_clean(n1,2:4)-node_positions_clean(n2,2:4)).^2,2));
    % for every element in edge list find d > dr
    r = find(edge_list(:,3) > dL);
    if ~isempty(r)
        n1 = edge_list(r,1);
        n2 = edge_list(r,2);
        d = edge_list(r,3);
        new_parts_num = [(tot_part_number + 1):(tot_part_number+length(r))]';
        tot_part_number = tot_part_number + length(r);
        edge_list(end+1:end+length(r),:) = [n1 new_parts_num d/2];
        edge_list(end+1:end+length(r),:) = [new_parts_num n2 d/2];
        pos_new = (node_positions(n1,2:4) + node_positions(n2,2:4))/2;
        node_positions(new_parts_num,:) = [new_parts_num pos_new];
        % if tot_part_number ~= length(
        edge_list(r,:) = [];
    end
end
end % end expand_domain

function [active_part_pos,active_part_angle,active_part_number,n_active_part, inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part] = check_reactivated_particles(active_part_pos,active_part_angle,active_part_number,n_active_part, inactive_part_pos,inactive_part_angle,inactive_part_number,n_inactive_part,annihil_radius,node_positions,sensing_angle)
if n_inactive_part > 0
    % trace_pos = node_positions(:,2:end);
    %
    all_pos = node_positions(:,2:end);
    for i = n_inactive_part:-1:1
        in_part_pos = inactive_part_pos(i,2:4);
        in_part_angle = inactive_part_angle(i,2:end);

        inds = (abs(all_pos(:,1)-in_part_pos(1)) <= annihil_radius & abs(all_pos(:,2)-in_part_pos(2)) <= annihil_radius & abs(all_pos(:,3)-in_part_pos(3)) <= annihil_radius);
        trace_pos = all_pos(inds,:);
        %         OPTIM
        dists = pt_to_pts_dist(in_part_pos,trace_pos);
        %         NON OPTIM
        %         dists = pdist2(in_part_pos,trace_pos);

        v = trace_pos - in_part_pos;

        w = [sin(in_part_angle(1))*cos(in_part_angle(2)) sin(in_part_angle(1))*sin(in_part_angle(2)) cos(in_part_angle(1))];
        
        w = repmat(w,size(v,1),1);
        angles = atan2(vecnorm(cross(v,w,2),2,2),dot(v,w,2))';
        %         angles = calcAngleBetweenVectors(v, w);

        % dists > 0 is to ignore in_part_pos, which has angle == 0
        h = any(dists > 0 & dists < annihil_radius & abs(angles) < sensing_angle, 'all');

        if ~h
            %             disp('reactivated!')
            active_part_pos(end+1,:) = inactive_part_pos(i,:);
            active_part_angle(end+1,:) = inactive_part_angle(i,:);
            active_part_number(end+1) = inactive_part_number(i);
            n_active_part = n_active_part + 1;

            inactive_part_pos(i,:) = [];
            inactive_part_angle(i,:) = [];
            inactive_part_number(i) = [];
            n_inactive_part = n_inactive_part - 1;
        end
    end
end
end

function distance = pt_to_pts_dist(pt,list)
distance = hypot(pt(:, 1) - list(:, 1).', pt(:, 2) - list(:, 2).');  %distance between all points
if size(pt,1) > 1
    distance(logical(tril(ones(size(distance))))) = Inf;  %point below diagonal are symmetric of upper triangle. Also remove diagonal from minimum search
end
% [mindistance, location] = min(distance(:));
% [point1, point2] = ind2sub(size(distance), location);
end

function [R, Rm] = rotate_3D(V, mode, theta, u, angle_unit)
% rotate_3D : function to compute the rotation of a vector or an array of vectors in 2D or 3D space.
% Source:
% Nicosahedron (2022). Any 3D rotation (https://github.com/NicolasDouillet/rotate_3D/releases/tag/v2.6), GitHub. Retrieved April 12, 2022.
%
% Syntax
% R = rotate_3D(V, mode, theta);
% R = rotate_3D(V, mode, theta, u);
% R = rotate_3D(V, mode, theta, u, angle_unit);
% [R,Rm] = rotate_3D(V, mode, theta, u, angle_unit);
%
%
% Description
% R = rotate_3D(V, mode, theta) computes the vector R, which results
% from the rotation of V vector around one of the the basis vectors, which
% is choosen in the mode : 'x', 'y', or 'z'.
%
% R = rotate_3D(V, mode, theta, u) computes the vector R, which results
% from the rotation of V vector around u vector and of theta angle in radian.
%
% R = rotate_3D(V, mode, theta, u, angle_unit) uses angle_unit for theta
% unit (radian or degree).
%
% [R,Rm] = rotate_3D(V, mode, theta, u, angle_unit) also returns the
% rotation matrix.
%
% Important NB : in 2D -(xOy) plan- mandatory rotation axis is 'z'. It will
% be set as so by default if input is different. Also in 2D, in case u is missing it
% is automatically set to the origin [0,0]' by default.
%
% Input parsing
Ndim = size(V,1);
assert(nargin > 2, 'Not enough input arguments.');
assert(nargin < 6, 'Too many input arguments.');
assert(Ndim > 1 && Ndim < 4, 'Input argument V must have between one and three rows : 1 < size(V,1) <= 3.');
assert(strcmpi(mode,'x') || strcmpi(mode,'y') || strcmpi(mode,'z') || strcmpi(mode,'any'),...
       'Bad mode argument : mode must be a string in the set {''x'',''X'',''y'',''Y'',''z'',''Z'',''any'',''ANY''}.');
if nargin < 5    
    angle_unit = 'radian';    
    if nargin < 4        
        if Ndim == 2           
            u = [0,0]';    
        elseif Ndim == 3    
            switch mode   
                case {'x', 'X'}
                    u = [1 0 0]';  
                case {'y', 'Y'}
                    u = [0 1 0]';
                case {'z', 'Z'}
                    u = [0 0 1]';     
            end 
        end
    else
        assert(Ndim < 3 || ~strcmpi(mode,'any') || norm(u) > 0,'3D rotation axis u must not equal null vector.');  
    end
else
    assert(strcmpi(angle_unit,'radian') || strcmpi(angle_unit,'degree'),'angle_unit value must be either ''radian'' or ''degree''.');
    if strcmpi(angle_unit,'degree')
        theta = pi * theta / 180;
    end
end
% Body
% Rotation matrix construction and resulting rotated vector computation
switch Ndim
    case 2 % rotate around a point (2D vector) in (xOy) plan -> mandatory rotation axis is 'z' 
        Rm = [cos(theta) -sin(theta);
              sin(theta)  cos(theta)];
        W = V - u;
        R = Rm * W;
        R = R + u;                
    case 3
        switch mode
            case {'x', 'X'} % X axis rotation matrix ; u = i = [1 0 0]'
                Rm = [1          0           0;
                      0 cos(theta) -sin(theta);
                      0 sin(theta)  cos(theta)];
            case {'y', 'Y'} % Y axis rotation matrix ; u = j = [0 1 0]'
                Rm = [cos(theta)   0  -sin(theta);
                      0            1           0;
                      sin(theta)  0  cos(theta)];
            case {'z', 'Z'} % Z axis rotation matrix ; u = k = [0 0 1]'
                Rm = [cos(theta) -sin(theta) 0;
                      sin(theta)  cos(theta) 0;
                      0           0          1];
            case {'any', 'ANY'} % Any u axis rotation matrix
                u = u/norm(u);
                Rm = [u(1,1)^2+cos(theta)*(1-u(1,1)^2) (1-cos(theta))*u(1,1)*u(2,1)-u(3,1)*sin(theta) (1-cos(theta))*u(1,1)*u(3,1)+u(2,1)*sin(theta);
                      (1-cos(theta))*u(1,1)*u(2,1)+u(3,1)*sin(theta) u(2,1)^2+cos(theta)*(1-u(2,1)^2) (1-cos(theta))*u(2,1)*u(3,1)-u(1,1)*sin(theta);
                      (1-cos(theta))*u(1,1)*u(3,1)-u(2,1)*sin(theta) (1-cos(theta))*u(2,1)*u(3,1)+u(1,1)*sin(theta) u(3,1)^2+cos(theta)*(1-u(3,1)^2)];
            otherwise
                error('Bad mode argument : mode must be a string in the set {''x'',''X'',''y'',''Y'',''z'',''Z'',''any'',''ANY''}.');   
        end
        R = Rm * V;            
end
end % rotate_3D

function y  = laprnd(m, n, mu, sigma)
%LAPRND generate i.i.d. random number drawn from Leplace distribution
%   with mean mu and standard deviation sigma. 
%
% Source: 
% Elvis Chen (2022). Laplacian random number generator (https://www.mathworks.com/matlabcentral/fileexchange/13705-laplacian-random-number-generator), MATLAB Central File Exchange. Retrieved August 24, 2022.
% 
% Syntax
%   mu      : mean
%   sigma   : standard deviation
%   [m, n]  : the dimension of y.
%   Default mu = 0, sigma = 1. 
%   For more information, refer to
%   http://en.wikipedia.org./wiki/Laplace_distribution

%   Author  : Elvis Chen (bee33@sjtu.edu.cn)
%   Date    : 01/19/07
%Check inputs
if nargin < 2
    error('At least two inputs are required');
end

if nargin == 2
    mu = 0; sigma = 1;
end

if nargin == 3
    sigma = 1;
end

% Generate Laplacian noise
u = rand(m, n)-0.5;
b = sigma / sqrt(2);
y = mu - b * sign(u).* log(1- 2* abs(u));
end % laprnd