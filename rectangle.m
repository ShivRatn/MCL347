
function fem_heat_conduction()
    % Parameters
    W = 1.0;  % Width of the domain
    H = 1.0;  % Height of the domain
    T1 = 100; % Base temperature
    
    % Mesh parameters
    nx = 40;  % Number of elements in x-direction
    ny = 40;  % Number of elements in y-direction
    
    % Generate mesh
    [p, e, t, dx, dy] = generate_mesh(W, H, nx, ny);
    
    % Number of nodes
    n_nodes = size(p, 2);
    
    % Assemble stiffness matrix and load vector
    [K, F] = assemble_system(p, t, n_nodes);
    
    % Apply boundary conditions
    [K, F] = apply_boundary_conditions(K, F, p, W, H, T1);
    
    % Solve the system
    T = K \ F;
    
    % Plotting
    plot_solution(p, t, T, W, H, T1);
    analytical_sol(W, H, T1, nx, ny);

end

function [p, e, t, dx, dy] = generate_mesh(W, H, nx, ny)
    % Generate a rectangular mesh
    dx = W / nx;
    dy = H / ny;
    
    % Create node coordinates
    [X, Y] = meshgrid(0:dx:W, 0:dy:H);
    p = [X(:)'; Y(:)'];
    
    % Create triangular elements
    t = [];
    for j = 1:ny
        for i = 1:nx
            % Node indices
            n1 = i + (j-1)*(nx+1);
            n2 = i+1 + (j-1)*(nx+1);
            n3 = i + j*(nx+1);
            n4 = i+1 + j*(nx+1);
            
            % Add two triangles for each rectangle
            t = [t, [n1; n3; n2]];  % Lower triangle
            t = [t, [n2; n3; n4]];  % Upper triangle
        end
    end
    e = [];

    t = t(1:3, :);
end

function [K, F] = assemble_system(p, t, n_nodes)
    % Initialize stiffness matrix and load vector
    K = sparse(n_nodes, n_nodes);
    F = zeros(n_nodes, 1);
    
    % Number of elements
    n_el = size(t, 2);
    
    for el = 1:n_el
        % Get node indices for this element
        nodes = t(1:3, el);
        
        % Get node coordinates
        x_coords = p(1, nodes);
        y_coords = p(2, nodes);
        
        % Calculate element matrices
        [K_el, F_el] = element_matrices(x_coords, y_coords);
        
        % Assemble into global system
        for i = 1:3
            for j = 1:3
                K(nodes(i), nodes(j)) = K(nodes(i), nodes(j)) + K_el(i, j);
            end
            F(nodes(i)) = F(nodes(i)) + F_el(i);
        end
    end
end

function [K_el, F_el] = element_matrices(x, y)

    
    % Area of the element
    area = 0.5 * abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)));
    
    % Derivatives of shape functions
    b = [y(2)-y(3), y(3)-y(1), y(1)-y(2)];
    c = [x(3)-x(2), x(1)-x(3), x(2)-x(1)];
    
    % Element stiffness matrix
    K_el = zeros(3, 3);
    for i = 1:3
        for j = 1:3
            K_el(i,j) = (b(i)*b(j) + c(i)*c(j)) / (4*area);
        end
    end
    
    % Element force vector (zero for Laplace equation)
    F_el = zeros(3, 1);
end

function [K, F] = apply_boundary_conditions(K, F, p, W, H, T1)
    % Apply Dirichlet boundary conditions
    n_nodes = size(p, 2);
    tol = 1e-10;
    
    for i = 1:n_nodes
        x = p(1, i);
        y = p(2, i);
        
        % Left boundary
        if abs(x) < tol
            K(i,:) = 0;
            K(i,i) = 1;
            F(i) = T1;
        end
        
        % Right boundary
        if abs(x - W) < tol
            K(i,:) = 0;
            K(i,i) = 1;
            F(i) = T1;
        end
        
        % Top boundary (sinusoidal)
        if abs(y - H) < tol
            K(i,:) = 0;
            K(i,i) = 1;
            F(i) = T1 * (1 + sin(pi * x / W));
        end

    end
end

function plot_solution(p, t, T, W, H, T1)
    % Create figure
    figure('Position', [100, 100, 800, 700]);
    
    % Create a triangulation object for interpolation
    tri = triangulation(t(1:3,:)', p(1,:)', p(2,:)');
    
    % Calculate temperature gradients (heat flux)
    [px, py] = calculateGradients(tri, T);
    
    % Create a grid for visualization
    n_grid = 100;
    [X, Y, Z] = createGrid(tri, T, n_grid);
    
    % Calculate gradients on the grid
    [X_grad, Y_grad, U, V] = createGradientGrid(tri, px, py, 20);
    
    % Plot the contours with gradient vectors
    [C, h] = contourf(X, Y, Z, 15);
    hold on;
    
    % Add contour lines
    contour(X, Y, Z, 15, 'LineWidth', 1, 'LineColor', 'k');
    

;
    
    % Adjust plot settings
    xlabel('x-coordinate', 'FontSize', 12);
    ylabel('y-coordinate', 'FontSize', 12);
    title('Temperature Field with Heat Flux Vectors', 'FontSize', 14);
    colorbar;
    colormap('jet');
    axis equal tight;
    grid on;
    

    % Add legend
    h_quiver = quiver(0.8*W, 0.1*H, -0.1, 0, 'w', 'LineWidth', 1.2);
   
end

function [px, py] = calculateGradients(tri, T)
    % Calculate the temperature gradient at each node
    
    % Get the triangulation points and connectivity
    P = tri.Points;
    t = tri.ConnectivityList;
    n_nodes = size(P, 1);
    n_elements = size(t, 1);
    
    % Initialize arrays to store gradients
    px = zeros(n_nodes, 1);
    py = zeros(n_nodes, 1);
    areas = zeros(n_nodes, 1);
    
    % Calculate element-wise gradients and project to nodes
    for e = 1:n_elements
        % Get the nodes for this element
        nodes = t(e, :);
        
        % Get coordinates
        x = P(nodes, 1);
        y = P(nodes, 2);
        
        % Get temperatures at these nodes
        T_el = T(nodes);
        
        % Calculate area of the element
        area = 0.5 * abs((x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)));
        
        % Calculate basis function derivatives
        b = [y(2)-y(3), y(3)-y(1), y(1)-y(2)];
        c = [x(3)-x(2), x(1)-x(3), x(2)-x(1)];
        
        % Calculate temperature gradients in this element
        dTdx = b * T_el / (2 * area);
        dTdy = c * T_el / (2 * area);
        
        % Accumulate gradients to nodes (weighted by element area)
        for i = 1:3
            px(nodes(i)) = px(nodes(i)) + dTdx * area;
            py(nodes(i)) = py(nodes(i)) + dTdy * area;
            areas(nodes(i)) = areas(nodes(i)) + area;
        end
    end
    
    % Normalize by total area
    for i = 1:n_nodes
        if areas(i) > 0
            px(i) = px(i) / areas(i);
            py(i) = py(i) / areas(i);
        end
    end
end

function [X, Y, T] = analytical_sol(W, H, T1, nx, ny)


    % Create mesh
    x = linspace(0, W, nx);
    y = linspace(0, H, ny);
    [X, Y] = meshgrid(x, y);

    % Compute analytical solution
    T = T1 + T1 * sin(pi * X / W) .* cosh(pi * Y / W) / cosh(pi * H / W);
    figure;
    contourf(X, Y, T, 20);
    colorbar;
    colormap('parula');
    title('Analytical Solution: Long Rectangular Bar');
    xlabel('x'); ylabel('y');
    grid on;
end

function [X, Y, Z] = createGrid(tri, T, n)
    % Create a regular grid for visualization
    
    % Get the bounding box
    P = tri.Points;
    xmin = min(P(:,1));
    xmax = max(P(:,1));
    ymin = min(P(:,2));
    ymax = max(P(:,2));
    
    % Create a regular grid
    [X, Y] = meshgrid(linspace(xmin, xmax, n), linspace(ymin, ymax, n));
    
    % Interpolate temperature values at the grid points
    F = scatteredInterpolant(P(:,1), P(:,2), T);
    Z = F(X, Y);
end

function [X, Y, U, V] = createGradientGrid(tri, px, py, n)
    % Create a coarser grid for gradient visualization
    
    % Get the bounding box
    P = tri.Points;
    xmin = min(P(:,1));
    xmax = max(P(:,1));
    ymin = min(P(:,2));
    ymax = max(P(:,2));
    
    % Create a regular but coarser grid for the gradient vectors
    [X, Y] = meshgrid(linspace(xmin, xmax, n), linspace(ymin, ymax, n));
    
    % Interpolate gradient values at the grid points
    Fx = scatteredInterpolant(P(:,1), P(:,2), px);
    Fy = scatteredInterpolant(P(:,1), P(:,2), py);
    
    U = Fx(X, Y);
    V = Fy(X, Y);
    
    % Normalize vectors for better visualization
    mag = sqrt(U.^2 + V.^2);
    max_mag = max(mag(:));
    if max_mag > 0
        U = U / max_mag;
        V = V / max_mag;
    end
end