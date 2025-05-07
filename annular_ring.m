function axisymmetric_heat_transfer()
    % Parameters
    r1 = 0.5;       % Inner radius [m]
    r2 = 1.0;       % Outer radius [m]
    T1 = 100;       % Inner temperature [°C]
    amplitude = 50; % Angular variation amplification
    nr = 80;        % Radial elements
    ntheta = 120;   % Angular elements

    % Generate structured polar mesh
    [nodes, elements, inner_nodes, outer_nodes] = generate_mesh(r1, r2, nr, ntheta);
    
    % Assemble axisymmetric stiffness matrix
    K = assemble_axisymmetric_stiffness(nodes, elements);
    
    % Apply boundary conditions
    [K, F] = apply_dirichlet_bcs(K, nodes, inner_nodes, outer_nodes, ntheta, T1, amplitude);
    
    % Solve system
    T = K \ F;
    
    % Post-processing
    verify_boundary_conditions(nodes, T, inner_nodes, outer_nodes, T1, amplitude);
    plot_comparison(nodes, T, r1, r2, amplitude);
    analytical_solution_plot_enhanced();
end

function [nodes, elements, inner_nodes, outer_nodes] = generate_mesh(r1, r2, nr, ntheta)
    % Create structured polar mesh with periodic boundary
    dr = (r2 - r1)/nr;
    dtheta = 2*pi/ntheta;
    
    % Node generation (polar -> Cartesian)
    nodes = [];
    for i = 0:nr
        r = r1 + i*dr;
        for j = 0:ntheta-1
            theta = j*dtheta;
            nodes(end+1,:) = [r*cos(theta), r*sin(theta)]; %#ok<AGROW>
        end
    end

    % Triangular element connectivity
    elements = [];
    for i = 0:nr-1
        for j = 0:ntheta-1
            j_next = mod(j+1, ntheta);
            n1 = i*ntheta + j + 1;
            n2 = i*ntheta + j_next + 1;
            n3 = (i+1)*ntheta + j + 1;
            n4 = (i+1)*ntheta + j_next + 1;
            
            elements(end+1,:) = [n1, n2, n4]; %#ok<AGROW>
            elements(end+1,:) = [n1, n4, n3]; %#ok<AGROW>
        end
    end

    % Identify boundary nodes
    inner_nodes = 1:ntheta;
    outer_nodes = (nr*ntheta + 1):(nr+1)*ntheta;
end

function K = assemble_axisymmetric_stiffness(nodes, elements)
    % Assemble stiffness matrix with axisymmetric terms
    n_nodes = size(nodes, 1);
    K = sparse(n_nodes, n_nodes);
    
    for el = 1:size(elements, 1)
        conn = elements(el,:);
        coords = nodes(conn,:);
        [theta, r] = cart2pol(coords(:,1), coords(:,2));
        
        % Axisymmetric element matrix
        Ke = axisymmetric_element(r, theta);
        
        % Assemble
        K(conn, conn) = K(conn, conn) + Ke;
    end
end

function Ke = axisymmetric_element(r, theta)
    % Axisymmetric triangular element formulation
    x = r.*cos(theta);
    y = r.*sin(theta);
    
    % Element area
    A = 0.5*abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)));
    r_avg = mean(r);
    
    % Standard stiffness matrix
    b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
    c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];
    Ke = (b*b' + c*c')/(4*A);
    
    % Axisymmetric scaling (2πr term)
    Ke = 2*pi*r_avg * Ke;
end

function [K, F] = apply_dirichlet_bcs(K, nodes, inner_nodes, outer_nodes, ntheta, T1, amplitude)
    % Apply essential boundary conditions
    n_nodes = size(K,1);
    F = zeros(n_nodes,1);
    
    % Dirichlet values with exact angular alignment
    dtheta = 2*pi/ntheta;
    theta = (0:ntheta-1)*dtheta;
    T_outer = 200 + amplitude*sin(theta);
    
    % Create boundary value vector
    T_vals = zeros(n_nodes,1);
    T_vals(inner_nodes) = T1;
    T_vals(outer_nodes) = T_outer;
    
    % Adjust force vector
    F = F - K*T_vals;
    
    % Apply boundary conditions
    K(outer_nodes,:) = 0; K(:,outer_nodes) = 0;
    K(inner_nodes,:) = 0; K(:,inner_nodes) = 0;
    K(inner_nodes, inner_nodes) = speye(length(inner_nodes));
    K(outer_nodes, outer_nodes) = speye(length(outer_nodes));
    
    F(inner_nodes) = T1;
    F(outer_nodes) = T_outer;
end

function plot_comparison(nodes, T, r1, r2, amplitude)
    % Create analytical solution
    [Xa, Ya, Ta] = analytical_solution(r1, r2, amplitude);
    
    % Interpolate FEM solution
    Ti = griddata(nodes(:,1), nodes(:,2), T, Xa, Ya, 'natural');
    
    % Plot comparison
    figure('Position', [100, 100, 1200, 500])
    
    % subplot(1,2,1)
    contourf(Xa, Ya, Ti, 40, 'EdgeColor', 'none')
    title('Numerical Solution Using FEM')
    caxis([100, 250])
    

    
    colormap(jet)
    colorbar
    axis equal tight
end

function [X, Y, T] = analytical_solution(r1, r2, amplitude)
    % Analytical solution components
    nr = 100;
    ntheta = 100;
    r = linspace(r1, r2, nr);
    theta = linspace(0, 2*pi, ntheta);
    [R, THETA] = meshgrid(r, theta);
    
    % Axisymmetric base solution
    T0 = ((200 - 100)/log(r2/r1)) * log(R/r1) + 100;
    
    % Angular variation component
    T1 = amplitude * (r2/(r2^2 - r1^2)) * (R - (r1^2)./R) .* sin(THETA);
    
    X = R.*cos(THETA);
    Y = R.*sin(THETA);
    T = T0 + T1;
end

function verify_boundary_conditions(nodes, T, inner_nodes, outer_nodes, T1, amplitude)
    % Verify inner BC
    inner_error = max(abs(T(inner_nodes) - T1));
    fprintf('Inner boundary error: %.2e°C\n', inner_error)
    
    % Verify outer BC
    [theta, ~] = cart2pol(nodes(outer_nodes,1), nodes(outer_nodes,2));
    expected = 200 + amplitude*sin(theta);
    outer_error = max(abs(T(outer_nodes) - expected));
    fprintf('Outer boundary error: %.2e°C\n', outer_error)
    
    % Plot BC verification
    figure
    subplot(1,2,1)
    polarplot(theta, T(outer_nodes), 'ro', 'MarkerFaceColor', 'r')
    hold on
    polarplot(theta, expected, 'k--')
    title('Outer Boundary Condition Verification')
    
    subplot(1,2,2)
    plot(theta, T(outer_nodes), 'ro', 'DisplayName', 'FEM')
    hold on
    plot(theta, expected, 'k-', 'DisplayName', 'Analytical')
    xlabel('\theta (rad)')
    ylabel('Temperature (°C)')
    legend()
    grid on
end

function analytical_solution_plot_enhanced()
    % Parameters
    r1 = 0.5;   % Inner radius
    r2 = 1.0;   % Outer radius
    T1 = 100;   % Inner temperature (constant)
    amplitude = 50; % Amplify the angular variation (adjust this value)
    nr = 100;   % Radial resolution
    ntheta = 100; % Angular resolution

    % Generate grid
    r = linspace(r1, r2, nr);
    theta = linspace(0, 2*pi, ntheta);
    [R, THETA] = meshgrid(r, theta);

    % Compute analytical solution with amplified angular term
    T0 = ((200 - T1)/log(r2/r1)) * log(R/r1) + T1;
    T1_angular = amplitude * (r2/(r2^2 - r1^2)) * (R - (r1^2)./R) .* sin(THETA);
    T_total = T0 + T1_angular;

    % Plot
    figure;
    X = R .* cos(THETA);
    Y = R .* sin(THETA);
    contourf(X, Y, T_total, 30, 'LineColor', 'none');
    colorbar;
    colormap('jet');
    
    caxis([min(T_total(:)), max(T_total(:))]); 
    
    axis equal tight;
    title('Temperature Distribution with Angular Variation (Analytical solution)');
end