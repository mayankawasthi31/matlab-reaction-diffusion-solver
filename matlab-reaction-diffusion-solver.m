function varargout = read_inp_file_from_GMSD()
    % Domain limits
    x_start = 0; x_end = 3; y_start = 0; y_end = 2;

    % Open input file
    fid = fopen('BOOK.inp','r');
    if fid == -1, error('Could not open the inp file'); end

    % Skip header (assumed 9 lines)
    for i = 1:9, fgetl(fid); end

    % Read node data
    tot_dofs = 12;
    node_id = zeros(1, tot_dofs);
    x_glbl  = zeros(1, tot_dofs);
    y_glbl  = zeros(1, tot_dofs);
    for i = 1:tot_dofs
        tmp = sscanf(fgetl(fid), '%d,%e,%e', 3);
        node_id(i) = tmp(1);
        x_glbl(i)  = tmp(2);
        y_glbl(i)  = tmp(3);
    end
    fgetl(fid);  % skip blank line

    % Read element connectivity
    tot_elem = 6;
    elem_start = 1;
    ndes_per_elmn = 4;

    elem_connect = zeros(tot_elem, ndes_per_elmn);
    for ii = elem_start:tot_elem
        tmp1 = sscanf(fgetl(fid), '%d,%d,%d,%d,%d', 5);
        elem_connect(ii,:) = tmp1(2:5);
    end
    for i = 1:18, fgetl(fid); end  % skip lines

    % Initial & Boundary definitions
    tot_surfaces = 4;
    bc_types = {'N','D','D','N'};
    flux_values = {{'0','-30'},{'0','0'},{'0','0'},{'-20','0'}};
    dirichlet_values = {'0','0','cos(pi*x/6)','0'};
    conc_old = zeros(tot_dofs,1);

    % Time-step & Theta Scheme
    theta = 0.5; deltat = 0.1; tn = 10;
    t_plot = min(5, tn);

    % Parse boundary info and compute normals
    Boundary_info = cell(1, tot_surfaces);
    Normals_info  = cell(1, tot_surfaces);
    for s = 1:tot_surfaces
        % Surface nodes line
        clean = regexprep(fgetl(fid),'[,]+',' ');
        vals  = sscanf(clean,'%d');
        surf_nodes = vals(1):vals(3):vals(2);
        fgetl(fid);
        % Surface elements line
        clean = regexprep(fgetl(fid),'[,]+',' ');
        vals  = sscanf(clean,'%d');
        if numel(vals)==2, vals(3)=1; end
        surf_elems = vals(1):vals(3):vals(2);
        fgetl(fid);

        % Preallocate to avoid resizing
        maxEdges = numel(surf_elems) * 4;
        b_info = zeros(maxEdges, 3);
        cnt = 0;
        for e = surf_elems
            con   = elem_connect(e,:);
            edges = [con(1) con(2); con(2) con(3); con(3) con(4); con(4) con(1)];
            for j = 1:4
                if all(ismember(edges(j,:), surf_nodes))
                    cnt = cnt + 1;
                    b_info(cnt, :) = [e, edges(j,:)];
                end
            end
        end
        if cnt == 0
            b_info = [];
        else
            b_info = b_info(1:cnt, :);
        end
        Boundary_info{s} = b_info;

        % Compute normals
        n_info = zeros(size(b_info,1), 2);
        for k = 1:size(b_info,1)
            e   = b_info(k,1);
            n1  = b_info(k,2);
            n2  = b_info(k,3);
            x1 = x_glbl(n1); y1 = y_glbl(n1);
            x2 = x_glbl(n2); y2 = y_glbl(n2);
            xc = mean(x_glbl(elem_connect(e,:)));
            yc = mean(y_glbl(elem_connect(e,:)));
            cand1 = [y2 - y1, -(x2 - x1)];
            cand2 = -cand1;
            mid   = 0.5 * [x1 + x2, y1 + y2];
            if dot(cand1, [xc, yc] - mid) < 0
                nvec = cand1;
            else
                nvec = cand2;
            end
            n_info(k,:) = nvec / norm(nvec);
        end
        Normals_info{s} = n_info;
    end
    fclose(fid);

    %% Gaussian Quadrature Setup
    NGPS = 4;
    tgp  = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3), -1/sqrt(3)];
    sgp  = [-1/sqrt(3), -1/sqrt(3), 1/sqrt(3), 1/sqrt(3)];
    wts  = [1, 1, 1, 1];

    total_entries = tot_elem * ndes_per_elmn^2;
    rows_K = zeros(total_entries,1);
    cols_K = zeros(total_entries,1);
    vals_K = zeros(total_entries,1);
    rows_C = zeros(total_entries,1);
    cols_C = zeros(total_entries,1);
    vals_C = zeros(total_entries,1);
    entry   = 1;
    q_global = zeros(tot_dofs,1);

    %% Assemble global matrices
    for elem = 1:tot_elem
        con  = elem_connect(elem,:);
        xnod = x_glbl(con)'; ynod = y_glbl(con)';
        Ke   = zeros(4); Ce   = zeros(4);
        D_vec = 15 + 0*xnod + 0*ynod;
        K_vec = 0  + 0*xnod + 0*ynod;
        for gp = 1:NGPS
            t = tgp(gp); s = sgp(gp); w = wts(gp);
            N    = 1/4*[(1-s)*(1-t); (1-s)*(1+t); (1+s)*(1+t); (1+s)*(1-t)];
            dNdt = 1/4*[-(1-s); (1-s); (1+s); -(1+s)];
            dNds = 1/4*[-(1-t); -(1+t); (1+t);  (1-t)];
            dxdt = dNdt' * xnod; dydt = dNdt' * ynod;
            dxds = dNds' * xnod; dyds = dNds' * ynod;
            J    = [dxdt, dydt; dxds, dyds]; invJ = inv(J);
            D_bar = N' * D_vec;
            K_bar = N' * K_vec;
            Bbar  = [dNdt'; dNds'; N'];
            Abr   = diag([D_bar, D_bar, K_bar]);
            Ke = Ke + (Bbar' * [invJ, [0;0]; 0,0,1]' * Abr * [invJ, [0;0]; 0,0,1] * Bbar) * (det(J)*w);
            Ce = Ce + (N * N') * (det(J)*w);
        end
        for i = 1:4
            for j = 1:4
                rows_K(entry) = con(i); cols_K(entry) = con(j); vals_K(entry) = Ke(i,j);
                rows_C(entry) = con(i); cols_C(entry) = con(j); vals_C(entry) = Ce(i,j);
                entry = entry + 1;
            end
        end
    end
    K_global = sparse(rows_K, cols_K, vals_K, tot_dofs, tot_dofs);
    C_global = sparse(rows_C, cols_C, vals_C, tot_dofs, tot_dofs);

    %% Boundary integration for q_global
    ngps_line   = 2;
    gauss_pts_l = [-1/sqrt(3), 1/sqrt(3)];
    gauss_wts_l = [1, 1];
    for s_idx = 1:tot_surfaces
        bc_type = bc_types{s_idx};
        flux_ex = flux_values{s_idx};
        b_info  = Boundary_info{s_idx};
        n_info  = Normals_info{s_idx};
        for i_edge = 1:size(b_info,1)
            elem  = b_info(i_edge,1);
            nodeA = b_info(i_edge,2);
            nodeB = b_info(i_edge,3);
            normal= n_info(i_edge,:);
            con   = elem_connect(elem,:);
            xnod  = x_glbl(con)'; ynod = y_glbl(con)';
            posA  = find(con==nodeA);
            posB  = find(con==nodeB);
            sortedPos = sort([posA, posB]);
            if isequal(sortedPos,[1 2])
                fixed = -1; freeVar = 't';
            elseif isequal(sortedPos,[2 3])
                fixed =  1; freeVar = 's';
            elseif isequal(sortedPos,[3 4])
                fixed =  1; freeVar = 't';
            else
                fixed = -1; freeVar = 's';
            end
            Qe = zeros(4,1);
            if strcmp(bc_type,'N')
                for ip = 1:ngps_line
                    xi = gauss_pts_l(ip);
                    if freeVar=='t'
                        t_pt = xi; s_pt = fixed;
                    else
                        t_pt = fixed; s_pt = xi;
                    end
                    Nedge = zeros(4,1);
                    for kN = 1:4
                        ti = [-1,1,1,-1]'; si = [-1,-1,1,1]';
                        Nedge(kN) = 1/4*(1 + ti(kN)*t_pt)*(1 + si(kN)*s_pt);
                    end
                    dNdt = 1/4*[-(1-s_pt); (1-s_pt); (1+s_pt); -(1+s_pt)];
                    dNds = 1/4*[-(1-t_pt); -(1+t_pt); (1+t_pt);  (1-t_pt)];
                    dxdt = dNdt' * xnod; dydt = dNdt' * ynod;
                    dxds = dNds' * xnod; dyds = dNds' * ynod;
                    J = [dxdt, dydt; dxds, dyds];
                    scaling = (freeVar=='t') * norm(J(1,:)) + (freeVar=='s') * norm(J(2,:));
                    xm = Nedge' * xnod; ym = Nedge' * ynod;
                    x = xm; y = ym; % for eval
                    flux_x = eval(flux_ex{1});
                    flux_y = eval(flux_ex{2});
                    D_val  = 15 + 0*x + 0*y;
                    integrand = Nedge * (normal(1)*D_val*flux_x + normal(2)*D_val*flux_y);
                    Qe = Qe + integrand * scaling * gauss_wts_l(ip);
                end
            end
            for ii = 1:4
                q_global(con(ii)) = q_global(con(ii)) + Qe(ii);
            end
        end
    end

    %% Time-stepping solver:
    for s = 1:tot_surfaces
        if strcmp(bc_types{s},'D')
            expr = dirichlet_values{s};
            nodes= unique(Boundary_info{s}(:,2:3));
            for nn = nodes'
                x = x_glbl(nn);
                y = y_glbl(nn);
                conc_old(nn) = eval(expr);
            end
        end
    end
    K1 = C_global + theta * deltat * K_global;
    K2 = C_global - (1-theta) * deltat * K_global;
    Q1 = deltat * ((1-theta)*q_global + theta*q_global);
    nt = ceil(tn / deltat);
    conc_history = zeros(tot_dofs, nt+1);
    conc_history(:,1) = conc_old;
    for step = 1:nt
        conc_new = zeros(tot_dofs,1); known = false(tot_dofs,1);
        for s = 1:tot_surfaces
            if strcmp(bc_types{s},'D')
                expr  = dirichlet_values{s}; nodes = unique(Boundary_info{s}(:,2:3));
                for nn = nodes'
                    x = x_glbl(nn); y = y_glbl(nn);
                    conc_new(nn) = eval(expr); known(nn) = true;
                end
            end
        end
        rhs = K2 * conc_old + Q1;
        rhs(known) = conc_new(known);
        unknown = ~known;
        conc_new(unknown) = K1(unknown,unknown)\(rhs(unknown) - K1(unknown,known)*conc_new(known));
        conc_history(:,step+1) = conc_new; conc_old = conc_new;
    end

    %% Debug plot
    figure;
    patch('Faces', elem_connect, 'Vertices', [x_glbl', y_glbl'], ...
          'FaceVertexCData', conc_history(:,1), 'FaceColor','interp', 'EdgeColor','none');
    colorbar; title('Initial concentration at t=0');

    %% Animation
    figure;
    for step = 1:nt+1
        patch('Faces', elem_connect, 'Vertices', [x_glbl', y_glbl'], ...
              'FaceVertexCData', conc_history(:,step), 'FaceColor','interp', 'EdgeColor','none');
        colorbar; title(sprintf('Concentration at t=%.2f', (step-1)*deltat));
        drawnow;
    end

    %% Snapshot of t_plot 
    step_lo = floor(t_plot/deltat); step_hi = ceil(t_plot/deltat);
    alpha = (t_plot - step_lo*deltat)/deltat;
    conc_lo = conc_history(:, step_lo+1);
    conc_hi = conc_history(:, step_hi+1);
    conc_plot = (1-alpha)*conc_lo + alpha*conc_hi;
    figure;
    patch('Faces', elem_connect, 'Vertices', [x_glbl', y_glbl'], ...
          'FaceVertexCData', conc_plot, 'FaceColor','interp', 'EdgeColor','none');
    colorbar; title(sprintf('Concentration at t=%.2f s', t_plot));

    %% Assign to base
    assignin('base','K_global',K_global);
    assignin('base','C_global',C_global);
    assignin('base','q_global',q_global);
    assignin('base','conc_history',conc_history);

    if nargout>0
        varargout = {K_global,C_global,q_global,x_start,x_end,y_start,y_end,...
                    x_glbl,y_glbl,elem_connect,tot_elem,node_id,tot_dofs,elem_start,ndes_per_elmn};
    end
end