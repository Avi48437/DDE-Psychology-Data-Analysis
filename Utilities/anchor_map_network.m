function anchor_map_network(M)

% ANCHOR_MAP_NETWORK
%
% Compact horizontal hierarchical representation:
%
% A2 -> Positive / Neutral / Negative -> Domain -> A1 -> Anchor block
%
% Each A1 factor has one anchor block.
%
% Layer-2 interpretation:
%
%   Positive : B2 loading > 0
%   Neutral  : B2 loading = 0
%   Negative : B2 loading < 0
%
% The B2 intercept is intentionally not displayed.

if M.K2==0
    error('No active Layer-2 factors are available for the network plot.');
end

%% ============================================================
% Horizontal locations
%% ============================================================

x_A2 = 80;
x_sign = 260;
x_domain = 450;
x_A1 = 630;
x_anchor = 1060;

A2_width = 100;
sign_width = 125;
domain_width = 95;
A1_width = 90;
anchor_width = 690;

node_height = 42;

%% ============================================================
% Compact anchor blocks
%% ============================================================

anchor_line_height = 20;
anchor_padding = 10;

%% ============================================================
% Compact vertical spacing
%% ============================================================

factor_gap = 18;
domain_gap = 30;
sign_gap = 50;
A2_gap = 70;

top_margin = 30;
bottom_margin = 30;

%% ============================================================
% Domains
%% ============================================================

domain_order = unique(M.factor_domains,'stable');

n_domains = numel(domain_order);

domain_colors = lines(n_domains);

%% ============================================================
% Height of each A1 subtree
%% ============================================================

factor_height = zeros(M.K1,1);
anchor_block_height = zeros(M.K1,1);

for k = 1:M.K1

    n_anchor = ...
        sum(strlength(M.anchor_items(:,k))>0);

    if n_anchor==0
        n_anchor = 1;
    end

    anchor_block_height(k) = ...
        2*anchor_padding + ...
        n_anchor*anchor_line_height;

    factor_height(k) = ...
        max(node_height,anchor_block_height(k));

end

%% ============================================================
% Height required by each A2 hierarchy
%% ============================================================

A2_height = zeros(M.K2,1);

for h = 1:M.K2

    loading = M.B2(:,h+1);

    sign_index = { ...
        find(loading>0), ...
        find(loading==0), ...
        find(loading<0)};

    total_height = 0;
    n_sign_used = 0;

    for s = 1:3

        current_idx = sign_index{s};

        if isempty(current_idx)
            continue
        end

        sign_height = 0;
        n_domain_used = 0;

        for d = 1:n_domains

            factor_idx = current_idx( ...
                M.factor_domains(current_idx)==domain_order(d));

            if isempty(factor_idx)
                continue
            end

            domain_height = ...
                sum(factor_height(factor_idx));

            if numel(factor_idx)>1

                domain_height = ...
                    domain_height + ...
                    factor_gap*(numel(factor_idx)-1);

            end

            if n_domain_used>0
                sign_height = sign_height + domain_gap;
            end

            sign_height = ...
                sign_height + domain_height;

            n_domain_used = n_domain_used+1;

        end

        if n_sign_used>0
            total_height = total_height + sign_gap;
        end

        total_height = ...
            total_height + sign_height;

        n_sign_used = n_sign_used+1;

    end

    A2_height(h) = max(total_height,120);

end

%% ============================================================
% Total logical height
%% ============================================================

canvas_height = ...
    top_margin + ...
    sum(A2_height) + ...
    A2_gap*max(M.K2-1,0) + ...
    bottom_margin;

canvas_width = 1430;

%% ============================================================
% Create figure
%
% IMPORTANT:
% No Position / OuterPosition / WindowState is specified.
% MATLAB can therefore place the figure in the common
% Figure Container as its own tab.
%% ============================================================

fig = figure( ...
    'Name',char(M.dataset_name + " Anchor Map"), ...
    'Color','w', ...
    'NumberTitle','off');

ax = axes( ...
    'Parent',fig, ...
    'Units','normalized', ...
    'Position',[0.02 0.04 0.96 0.90]);

hold(ax,'on');

ax.XLim = [0 canvas_width];
ax.YLim = [0 canvas_height];

ax.YDir = 'reverse';
ax.Visible = 'off';
ax.Color = 'w';

%% ============================================================
% Title
%% ============================================================

title( ...
    ax, ...
    M.dataset_name + " Hierarchical Anchor Map", ...
    'FontSize',16, ...
    'FontWeight','bold', ...
    'Interpreter','none', ...
    'Color','k', ...
    'Visible','on');

%% ============================================================
% Sign definitions
%% ============================================================

sign_names = ["Positive","Neutral","Negative"];

sign_face = [ ...
    0.89 0.96 0.90
    0.94 0.94 0.94
    0.98 0.90 0.90];

sign_edge = [ ...
    0.20 0.55 0.25
    0.45 0.45 0.45
    0.75 0.20 0.20];

%% ============================================================
% Scale B2 edge widths
%% ============================================================

max_B2 = max(abs(M.B2(:,2:end)),[],'all');

if max_B2==0
    max_B2 = 1;
end

%% ============================================================
% Draw each A2 hierarchy
%% ============================================================

current_A2_top = top_margin;

for h = 1:M.K2

    loading = M.B2(:,h+1);

    sign_index = { ...
        find(loading>0), ...
        find(loading==0), ...
        find(loading<0)};

    A2_center_y = ...
        current_A2_top + ...
        A2_height(h)/2;

    %% --------------------------------------------------------
    % A2 node
    %% --------------------------------------------------------

    draw_box( ...
        ax, ...
        x_A2,A2_center_y, ...
        A2_width,node_height, ...
        "A2_" + string(h), ...
        [0.90 0.93 1.00], ...
        [0.20 0.30 0.65], ...
        11,true);

    current_sign_top = current_A2_top;
    n_sign_drawn = 0;

    %% ========================================================
    % Positive / Neutral / Negative
    %% ========================================================

    for s = 1:3

        current_idx = sign_index{s};

        if isempty(current_idx)
            continue
        end

        %% ----------------------------------------------------
        % Sign branch height
        %% ----------------------------------------------------

        sign_height = 0;
        n_domain_used = 0;

        for d = 1:n_domains

            factor_idx = current_idx( ...
                M.factor_domains(current_idx)==domain_order(d));

            if isempty(factor_idx)
                continue
            end

            domain_height = ...
                sum(factor_height(factor_idx));

            if numel(factor_idx)>1

                domain_height = ...
                    domain_height + ...
                    factor_gap*(numel(factor_idx)-1);

            end

            if n_domain_used>0
                sign_height = sign_height + domain_gap;
            end

            sign_height = ...
                sign_height + domain_height;

            n_domain_used = n_domain_used+1;

        end

        if n_sign_drawn>0
            current_sign_top = current_sign_top + sign_gap;
        end

        sign_center_y = ...
            current_sign_top + sign_height/2;

        %% ----------------------------------------------------
        % A2 -> sign
        %% ----------------------------------------------------

        draw_elbow( ...
            ax, ...
            x_A2+A2_width/2,A2_center_y, ...
            x_sign-sign_width/2,sign_center_y, ...
            sign_edge(s,:),1.9);

        %% ----------------------------------------------------
        % Sign box
        %% ----------------------------------------------------

        draw_box( ...
            ax, ...
            x_sign,sign_center_y, ...
            sign_width,node_height, ...
            sign_names(s), ...
            sign_face(s,:), ...
            sign_edge(s,:), ...
            10,true);

        %% ====================================================
        % Domains
        %% ====================================================

        current_domain_top = current_sign_top;
        n_domain_drawn = 0;

        for d = 1:n_domains

            factor_idx = current_idx( ...
                M.factor_domains(current_idx)==domain_order(d));

            if isempty(factor_idx)
                continue
            end

            domain_height = ...
                sum(factor_height(factor_idx));

            if numel(factor_idx)>1

                domain_height = ...
                    domain_height + ...
                    factor_gap*(numel(factor_idx)-1);

            end

            if n_domain_drawn>0

                current_domain_top = ...
                    current_domain_top + domain_gap;

            end

            domain_center_y = ...
                current_domain_top + domain_height/2;

            domain_color = domain_colors(d,:);

            domain_face = ...
                0.88 + 0.12*domain_color;

            %% ------------------------------------------------
            % Sign -> domain
            %% ------------------------------------------------

            draw_elbow( ...
                ax, ...
                x_sign+sign_width/2,sign_center_y, ...
                x_domain-domain_width/2,domain_center_y, ...
                sign_edge(s,:),1.4);

            %% ------------------------------------------------
            % Domain box
            %% ------------------------------------------------

            draw_box( ...
                ax, ...
                x_domain,domain_center_y, ...
                domain_width,node_height, ...
                domain_order(d), ...
                domain_face, ...
                domain_color, ...
                10,true);

            %% ================================================
            % Canonical A1 factors
            %% ================================================

            current_factor_top = current_domain_top;

            for ii = 1:numel(factor_idx)

                k = factor_idx(ii);

                factor_center_y = ...
                    current_factor_top + ...
                    factor_height(k)/2;

                %% --------------------------------------------
                % Domain -> A1
                %% --------------------------------------------

                if loading(k)==0

                    edge_color = sign_edge(2,:);
                    edge_width = 1.2;

                else

                    edge_color = sign_edge(s,:);

                    edge_width = ...
                        1.2 + ...
                        3*abs(loading(k))/max_B2;

                end

                draw_elbow( ...
                    ax, ...
                    x_domain+domain_width/2,domain_center_y, ...
                    x_A1-A1_width/2,factor_center_y, ...
                    edge_color,edge_width);

                %% --------------------------------------------
                % A1 box
                %% --------------------------------------------

                draw_box( ...
                    ax, ...
                    x_A1,factor_center_y, ...
                    A1_width,node_height, ...
                    M.factor_names(k), ...
                    domain_face, ...
                    domain_color, ...
                    10,true);

                %% ============================================
                % Build compact anchor block
                %% ============================================

                valid_anchor = find( ...
                    strlength(M.anchor_items(:,k))>0);

                anchor_lines = strings(0,1);

                for jj = 1:numel(valid_anchor)

                    j = valid_anchor(jj);

                    item_label = M.anchor_items(j,k);
                    question = M.anchor_questions(j,k);

                    if strlength(question)>0

                        anchor_lines(end+1,1) = ...
                            "• " + ...
                            question + ...
                            " [" + item_label + "]";

                    else

                        anchor_lines(end+1,1) = ...
                            "• " + item_label;

                    end

                end

                if isempty(anchor_lines)

                    anchor_text = ...
                        "• No positive anchor item";

                else

                    anchor_text = ...
                        strjoin(anchor_lines,newline);

                end

                %% --------------------------------------------
                % A1 -> anchor block
                %% --------------------------------------------

                draw_elbow( ...
                    ax, ...
                    x_A1+A1_width/2,factor_center_y, ...
                    x_anchor-anchor_width/2,factor_center_y, ...
                    domain_color,1.2);

                %% --------------------------------------------
                % Anchor block
                %% --------------------------------------------

                draw_anchor_block( ...
                    ax, ...
                    x_anchor,factor_center_y, ...
                    anchor_width,anchor_block_height(k), ...
                    anchor_text, ...
                    domain_face, ...
                    domain_color);

                %% --------------------------------------------
                % Next factor
                %% --------------------------------------------

                current_factor_top = ...
                    current_factor_top + ...
                    factor_height(k) + ...
                    factor_gap;

            end

            current_domain_top = ...
                current_domain_top + domain_height;

            n_domain_drawn = ...
                n_domain_drawn+1;

        end

        current_sign_top = ...
            current_sign_top + sign_height;

        n_sign_drawn = ...
            n_sign_drawn+1;

    end

    %% ========================================================
    % Separator between A2 factors
    %% ========================================================

    if h<M.K2

        separator_y = ...
            current_A2_top + ...
            A2_height(h) + ...
            A2_gap/2;

        plot( ...
            ax, ...
            [20 canvas_width-20], ...
            [separator_y separator_y], ...
            '--', ...
            'Color',[0.82 0.82 0.82], ...
            'LineWidth',0.8);

    end

    current_A2_top = ...
        current_A2_top + ...
        A2_height(h) + ...
        A2_gap;

end

%% ============================================================
% Force every text object to true black
%% ============================================================

set(findall(fig,'Type','text'),'Color','k');

end


%% ============================================================
% Draw centered node box
%% ============================================================

function draw_box(ax,x,y,w,h,label,face_color,edge_color,font_size,bold)

rectangle( ...
    ax, ...
    'Position',[x-w/2,y-h/2,w,h], ...
    'Curvature',0.10, ...
    'FaceColor',face_color, ...
    'EdgeColor',edge_color, ...
    'LineWidth',1.4);

if bold
    font_weight = 'bold';
else
    font_weight = 'normal';
end

text( ...
    ax, ...
    x,y, ...
    label, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontSize',font_size, ...
    'FontWeight',font_weight, ...
    'Interpreter','none', ...
    'Color','k');

end


%% ============================================================
% Draw compact anchor block
%% ============================================================

function draw_anchor_block(ax,x,y,w,h,label,face_color,edge_color)

rectangle( ...
    ax, ...
    'Position',[x-w/2,y-h/2,w,h], ...
    'Curvature',0.025, ...
    'FaceColor',face_color, ...
    'EdgeColor',edge_color, ...
    'LineWidth',1.2);

text( ...
    ax, ...
    x-w/2+14,y, ...
    label, ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','middle', ...
    'FontSize',9.5, ...
    'FontWeight','normal', ...
    'Interpreter','none', ...
    'Color','k');

end


%% ============================================================
% Draw horizontal elbow connector
%% ============================================================

function draw_elbow(ax,x1,y1,x2,y2,color,line_width)

x_mid = x1 + 0.45*(x2-x1);

plot( ...
    ax, ...
    [x1 x_mid x_mid x2], ...
    [y1 y1 y2 y2], ...
    '-', ...
    'Color',color, ...
    'LineWidth',line_width);

end