function fem_cylindrical_laplace()
    % Parameters
    R = 1.0;    % Radius of the cylinder
    L = 1.0;    % Half-length of the cylinder
    T1 = 200;   % Base temperature
    T0 = 100;   % Temperature amplitude
    m = 1;      % Angular mode number
    
    % Mesh parameters
    nr = 20;    % Number of elements in radial direction
    ntheta = 40; % Number of elements in angular direction
    nz = 20;    % Number of elements in axial direction
    
    % Generate mesh (2D cross-section in theta-z plane)
    [p, e, t] = generate_cylindrical_mesh(R, L, ntheta, nz);
    
    % Number of nodes
    n_nodes = size(p, 2);
    
    % Assemble stiffness matrix and load vector
    [K, F] = assemble_system(p, t, n_nodes, R);
    
    % Apply boundary conditions
    [K, F] = apply_boundary_conditions(K, F, p, R, L, T1, T0, m);
    
    % Solve the system
    T = K \ F;
    
    % Plotting
    plot_solution(p, t, T, R, L, T1, T0);
    plot_analytical_solution(R, L, T1, T0, m);
end

function [p, e, t] = generate_cylindrical_mesh(R, L, ntheta, nz)

    % Create node coordinates
    theta = linspace(0, 2*pi, ntheta+1);
    z = linspace(-L, L, nz+1);
    [Theta, Z] = meshgrid(theta(1:end-1), z); 
    
    % Convert to Cartesian coordinates for plotting (but solve in theta-z)
    p = [Theta(:)'; Z(:)'];
    
    % Create triangular elements
    t = [];
    for j = 1:nz
        for i = 1:ntheta
            % Node indices (periodic in theta)
            n1 = i + (j-1)*ntheta;
            n2 = mod(i, ntheta)+1 + (j-1)*ntheta;
            n3 = i + j*ntheta;
            n4 = mod(i, ntheta)+1 + j*ntheta;
            
            % Add two triangles for each rectangle
            t = [t, [n1; n3; n2]];  % Lower triangle
            t = [t, [n2; n3; n4]];  % Upper triangle
        end
    end
    
    % Create edge data (for boundary conditions)
    e = [];
end

function [K, F] = assemble_system(p, t, n_nodes, R)
    % Initialize stiffness matrix and load vector
    K = sparse(n_nodes, n_nodes);
    F = zeros(n_nodes, 1);
    
    % Number of elements
    n_el = size(t, 2);
    
    for el = 1:n_el
        % Get node indices for this element
        nodes = t(1:3, el);
        
        % Get node coordinates (theta, z)
        theta = p(1, nodes);
        z = p(2, nodes);
        
        % Calculate element matrices for cylindrical coordinates
        [K_el, F_el] = element_matrices_cylindrical(theta, z, R);
        
        % Assemble into global system
        for i = 1:3
            for j = 1:3
                K(nodes(i), nodes(j)) = K(nodes(i), nodes(j)) + K_el(i, j);
            end
            F(nodes(i)) = F(nodes(i)) + F_el(i);
        end
    end
end

function [K_el, F_el] = element_matrices_cylindrical(theta, z, R)
    
    % Area of the element in theta-z plane
    area = 0.5 * abs((theta(2)-theta(1))*(z(3)-z(1)) - (theta(3)-theta(1))*(z(2)-z(1)));
    
    % Derivatives of shape functions
    b = [z(2)-z(3), z(3)-z(1), z(1)-z(2)];
    c = [theta(3)-theta(2), theta(1)-theta(3), theta(2)-theta(1)]; 
    
    % Element stiffness matrix (cylindrical coordinates)
    K_el = zeros(3, 3);
    for i = 1:3
        for j = 1:3
            K_el(i,j) = (1/R^2 * b(i)*b(j) + c(i)*c(j)) / (4*area);
        end
    end
    
    % Element force vector (zero for Laplace equation)
    F_el = zeros(3, 1);
end

function [K, F] = apply_boundary_conditions(K, F, p, R, L, T1, T0, m)
    % Apply boundary conditions
    n_nodes = size(p, 2);
    tol = 1e-6;
    
    for i = 1:n_nodes
        theta = p(1, i);
        z = p(2, i);
        
        % Top boundary (z = L)
        if abs(z - L) < tol
            K(i,:) = 0;
            K(i,i) = 1;
            F(i) = T1 + T0*cos(m*theta);
        end
        
        % Bottom boundary (z = -L)
        if abs(z + L) < tol
            K(i,:) = 0;
            K(i,i) = 1;
            F(i) = T1 - T0*cos(m*theta);
        end
        % Enforce periodic BC in theta
    [K, F] = enforce_periodic_theta(K, F, p);

    end
end

function plot_solution(p, t, T, R, L, T1, T0)
    % Extract theta and z coordinates
    theta = p(1,:);
    z = p(2,:);
    
    % Convert to cylindrical surface coordinates
    x = R * cos(theta);
    y = R * sin(theta);

    figure('Name','3D Cylindrical Surface with Isotherms','Position',[100, 100, 1000, 800]);

    tri = triangulation(t(1:3,:)', x', y', z');
    h_surf = trisurf(tri, T, 'EdgeColor', 'none');
    colormap(jet(256));
    h_cb = colorbar;
    title(h_cb, 'Temperature');
    title(['Temperature Distribution (T_1 = ', num2str(T1), ', T_0 = ', num2str(T0), ')'], 'FontSize', 14);
    xlabel('r cos(\theta)', 'FontSize', 12);
    ylabel('r sin(\theta)', 'FontSize', 12);
    zlabel('z', 'FontSize', 12);
    view(3);
    axis equal tight;
    shading interp;
    hold on;

    % Interpolate T on a regular (theta,z) grid
    theta_grid = linspace(0, 2*pi, 400);  % finer for smooth contours
    z_grid = linspace(-L, L, 200);
    [TH, ZH] = meshgrid(theta_grid, z_grid);

    % Interpolate temperature values
    F_interp = scatteredInterpolant(theta', z', T, 'natural', 'none');
    T_grid = F_interp(TH, ZH);

    % Convert grid to (x,y) cylindrical surface
    X_grid = R * cos(TH);
    Y_grid = R * sin(TH);

    % Calculate isotherm levels based on temperature range and boundary conditions
    temp_min = min([min(T), T1-T0]);
    temp_max = max([max(T), T1+T0]);
    isotherm_levels = linspace(temp_min, temp_max, 15); % 15 isotherms

    [C3D, h3D] = contour3(X_grid, Y_grid, ZH, T_grid, isotherm_levels, 'k', 'LineWidth', 1.2);
    clabel(C3D, h3D, 'FontSize', 9, 'LabelSpacing', 500); % Label selected contours
    hold off;

    %--- 2D CONTOUR PLOT IN THETA-Z PLANE ---
    figure('Name','2D Contour in \theta-z Plane','Position',[150, 150, 900, 600]);
    
    [~, h_fill] = contourf(theta_grid, z_grid, T_grid, 40, 'LineColor', 'none');
    colormap(jet(256)); 
    h_cb = colorbar;
    title(h_cb, 'Temperature');
    hold on;
    
    % Add explicit isotherm lines on top of filled contours
    [C, h] = contour(theta_grid, z_grid, T_grid, isotherm_levels, 'k', 'LineWidth', 1.2);
    clabel(C, h, 'FontSize', 10, 'LabelSpacing', 500); % Add labels to contour lines
    
    xlabel('\theta', 'FontSize', 14);
    ylabel('z', 'FontSize', 14);
    title(['Isotherms in \theta-z Plane (T_1 = ', num2str(T1), ', T_0 = ', num2str(T0), ')'], 'FontSize', 16);
    set(gca, 'FontSize', 12);

    grid on;
    
    xticks([0 pi/2 pi 3*pi/2 2*pi]);
    xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    caxis([temp_min, temp_max]);
end


function [K, F] = enforce_periodic_theta(K, F, p)
    tol = 1e-8;
    theta_vals = p(1, :);
    z_vals = p(2, :);

    % Find nodes at theta ≈ 0 and theta ≈ 2π
    nodes_0 = find(abs(theta_vals - 0) < tol);
    nodes_2pi = find(abs(theta_vals - 2*pi) < tol | abs(theta_vals) > 2*pi - tol);

    % For each matching z-level, merge dofs
    for i = 1:length(nodes_0)
        n0 = nodes_0(i);
        z0 = z_vals(n0);

        % Find matching node at theta = 2π with same z
        match = find(abs(z_vals(nodes_2pi) - z0) < tol);
        if isempty(match)
            continue;
        end
        n2pi = nodes_2pi(match(1));  % Matching node at 2π

        % Merge n2pi into n0
        K(n0,:) = K(n0,:) + K(n2pi,:);
        F(n0) = F(n0) + F(n2pi);

        K(:,n0) = K(:,n0) + K(:,n2pi);

        % Zero out n2pi row and column
        K(n2pi,:) = 0;
        K(:,n2pi) = 0;
        K(n2pi,n2pi) = 1;
        F(n2pi) = F(n0); % Ensure continuity
    end
end

function plot_analytical_solution(R, L, T1, T0, m)
    
    % Default values if not provided
    if nargin < 1, R = 1.0; end
    if nargin < 2, L = 1.0; end
    if nargin < 3, T1 = 200; end
    if nargin < 4, T0 = 100; end
    if nargin < 5, m = 1; end
    
    % Create a fine grid for plotting
    theta_points = 100;
    z_points = 100;
    
    theta = linspace(0, 2*pi, theta_points);
    z = linspace(-L, L, z_points);
    [THETA, Z] = meshgrid(theta, z);
    
    % Calculate analytical solution on the grid
    hyperbolic_factor = sinh(m * Z / R) ./ sinh(m * L / R);
    T = T1 + T0 * cos(m * THETA) .* hyperbolic_factor;
    
    % Convert to cylindrical coordinates for 3D plotting
    X = R * cos(THETA);
    Y = R * sin(THETA);
    
    % Figure 1: 3D Cylindrical Surface Plot
    figure('Position', [100, 100, 1000, 800]);
    surf(X, Y, Z, T);
    colormap(jet(256));
    h_cb = colorbar;
    title(h_cb, 'Temperature');
    shading interp;
    axis equal;
    view(30, 20);
    title(sprintf('Analytical Solution: T(\\theta,z) = %.0f + %.0f cos(%d\\theta) \\cdot \\frac{sinh(%dz/%.1f)}{sinh(%dL/%.1f)}', ...
          T1, T0, m, m, R, m, R), 'FontSize', 14);
    xlabel('x = Rcos(\theta)', 'FontSize', 12);
    ylabel('y = Rsin(\theta)', 'FontSize', 12);
    zlabel('z', 'FontSize', 12);
    
    % Add isotherms (contour lines) on the 3D surface
    hold on;
    temp_min = min(T(:));
    temp_max = max(T(:));
    isotherm_levels = linspace(temp_min, temp_max, 15);
    [~, h_cont] = contour3(X, Y, Z, T, isotherm_levels, 'k', 'LineWidth', 1.2);
    hold off;
    
    % Figure 2: 2D Contour Plot in theta-z plane
    figure('Position', [150, 150, 900, 600]);
    
    % Create filled contour plot
    [~, h_fill] = contourf(theta, z, T, 40, 'LineColor', 'none');
    colormap(jet(256));
    h_cb = colorbar;
    title(h_cb, 'Temperature');
    hold on;
    
    % Add isotherms as explicit contour lines
    [C, h] = contour(theta, z, T, isotherm_levels, 'k', 'LineWidth', 1.2);
    clabel(C, h, 'FontSize', 10, 'LabelSpacing', 500);
    
    title(sprintf('Analytical Solution in \\theta-z Plane (T_1 = %.0f, T_0 = %.0f, m = %d)', ...
          T1, T0, m), 'FontSize', 16);
    xlabel('\theta', 'FontSize', 14);
    ylabel('z', 'FontSize', 14);
    grid on;
    
    % Special theta ticks at important angles
    xticks([0 pi/2 pi 3*pi/2 2*pi]);
    xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
    
end