classdef MarsBase
    methods(Static)
        function [N0,R0] = MakeInitialState(assumptions)
            
            % Generates initial binary populations and R0 resource conc
            % Outputs as struct but is converted to array later
            % Use struct2array() from essential_functions if needed

            M = sum(assumptions.MA); % total number of resources
            S_tot = sum(assumptions.SA)+assumptions.Sgen; % total number of consumers

            R0 = zeros(M,assumptions.n_wells);
            N0 = zeros(S_tot,assumptions.n_wells);

            % INDEXING

            % Creates species and family index
            F_index = strings(length(assumptions.SA)+1,1);
            S_index = struct();
            counter = 1;

            for i = 1:length(assumptions.SA)

                F_index(i) = 'F'+string(i);
                S_names = strings(assumptions.SA(i),1);

                for j = 1:assumptions.SA(i)
                    S_names(j) = 'S'+string(counter);
                    counter = counter+1;
                end

                S_index.('F'+string(i)) = S_names;
            end

            S_names = strings(assumptions.Sgen,1);
            for i = 1:assumptions.Sgen
                S_names(i) = 'S'+string(counter);
                counter = counter+1;
            end
            F_index(length(assumptions.SA)+1) = 'GEN';
            S_index.GEN = S_names;


            % Creates resource and type index

            T_index = strings(length(assumptions.MA),1);
            R_index = struct();
            counter = 1;

            for i = 1:length(assumptions.MA)

                T_index(i) = 'T'+string(i);
                R_names = strings(assumptions.MA(i),1);

                for j = 1:assumptions.MA(i)
                    R_names(j) = 'R'+string(counter);
                    counter = counter+1;
                end

                R_index.('T'+string(i)) = R_names;
            end

            % Creates well index
            well_index = strings(assumptions.n_wells,1);
            for i = 1:assumptions.n_wells
                well_index(i) = 'W'+string(i);
            end

            
            % COMPUTING VALUES

            % Sampling for N0, assigning given R0
            for i = 1:assumptions.n_wells
                N0(randsample(S_tot,assumptions.S),i) = 1;
                R0(:,i) = assumptions.R0;
            end

            % FORMATTING N0 AND R0

            % Formats each species family into table, searchable in N0_struct
            counter = 1;
            for i = 1:length(F_index)-1
                Fi = transpose(N0(counter:(counter+assumptions.SA(i)-1),:));
                N0_struct.('F'+string(i)) = array2table(Fi,'VariableNames',S_index.('F'+string(i)));
                counter = counter + assumptions.SA(i);
            end
            Gen = transpose(N0((S_tot-assumptions.Sgen+1):S_tot,:));
            N0_struct.('GEN') = array2table(Gen,'VariableNames',S_index.GEN);


            % Formats each resource type into table, searchable in R0_struct
            counter = 1;
            for i = 1:length(T_index)
                Ti = transpose(R0(counter:(counter+assumptions.MA(i)-1),:));
                R0_struct.('T'+string(i)) = array2table(Ti,'VariableNames',R_index.('T'+string(i)));
                counter = counter + assumptions.MA(i);
            end

            N0 = N0_struct;
            R0 = R0_struct;
        end

        function [c,D] = MakeMatrices(assumptions)
            
            % Generates consumption matrix c & metabolic byproducts matrix D
            % assumptions:  'sampling' - method for sampling c

            import essential_functions.drchrnd

            M = sum(assumptions.MA); % total number of resources
            S_tot = sum(assumptions.SA)+assumptions.Sgen; % total number of consumers

            c = []; % Will grow to (SxM)
            D = ones(M); % Will grow to (MxM)

            % Various sampling methods of c
            if strcmp(assumptions.sampling,'Gaussian')

                % Sampling for specialists per resource type
                for i = 1:length(assumptions.SA)
                    Fi = [];
                    for j = 1:length(assumptions.MA)

                        % Set distributions for specialists
                        % See Eq(4) in MiCRM paper
                        if i == j
                            c_mean = (assumptions.muc/M)*(1+assumptions.q*(M-assumptions.MA(j))/assumptions.MA(j));
                            c_var = (assumptions.sigc^2/M)*(1+assumptions.q*(M-assumptions.MA(j))/assumptions.MA(j));
                        else
                            c_mean = (assumptions.muc/M)*(1-assumptions.q);
                            c_var = (assumptions.sigc^2/M)*(1-assumptions.q);
                        end

                        FiTj = normrnd(c_mean,sqrt(c_var),assumptions.SA(i),assumptions.MA(j));
                        Fi = [Fi,FiTj];
                    end
                    c = [c;Fi];
                end

                if assumptions.Sgen > 0
                    % Set distributions for generalists
                    % See Eq(5) in MiCRM paper
                    c_mean = assumptions.muc/M;
                    c_var = assumptions.sigc^2/M;
                    Gen = [];

                    % Sampling for generalists per resource type
                    for j = 1:length(assumptions.MA)
                        GenTj = normrnd(c_mean,sqrt(c_var),assumptions.Sgen,assumptions.MA(j));
                        Gen = [Gen,GenTj];
                    end
                    c = [c;Gen];
                end

            elseif strcmp(assumptions.sampling,'Binary')
                assert(assumptions.muc < M*assumptions.c1, 'muc not attainable with given M and c1')

                % Sampling for specialists per resource type
                for i = 1:length(assumptions.SA)
                    Fi = [];
                    for j = 1:length(assumptions.MA)

                        % Set distributions for specialists
                        % See Eq(4) in MiCRM paper
                        if i == j
                            p = (assumptions.muc/(M*assumptions.c1))*(1+assumptions.q*(M-assumptions.MA(j))/assumptions.MA(j));
                        else
                            p = (assumptions.muc/(M*assumptions.c1))*(1-assumptions.q);
                        end

                        % Creates array of low consumption rate, generates binary
                        % array with probability p, 1s use high consumption rate
                        c0_array = assumptions.c0*ones(assumptions.SA(i),assumptions.MA(j));
                        unif = rand(assumptions.SA(i),assumptions.MA(j));
                        bin = unif<p;
                        FiTj = c0_array + (assumptions.c1-assumptions.c0)*bin;
                        Fi = [Fi,FiTj];
                    end
                    c = [c;Fi];
                end


                if assumptions.Sgen > 0
                    % Set distributions for generalists
                    % See Eq(5) in MiCRM paper
                    p = assumptions.muc/(M*assumptions.c1);
                    Gen = [];

                    % Sampling for generalists per resource type
                    for j = 1:length(assumptions.MA)

                        % Creates array of low consumption rate, generates binary
                        % array with probability p, 1s use high consumption rate
                        c0_array = assumptions.c0*ones(assumptions.Sgen,assumptions.MA(j));
                        unif = rand(assumptions.Sgen,assumptions.MA(j));
                        bin = unif<p;
                        GenTj = c0_array + (assumptions.c1-assumptions.c0)*bin;
                        Gen = [Gen,GenTj];
                    end
                    c = [c;Gen];
                end


            elseif strcmp(assumptions.sampling,'Gamma')

                % Sampling for specialists per resource type
                for i = 1:length(assumptions.SA)
                    Fi = [];
                    for j = 1:length(assumptions.MA)

                        % Set distributions for specialists
                        % See Eq(4) in MiCRM paper
                        if i == j
                            c_mean = (assumptions.muc/M)*(1+assumptions.q*(M-assumptions.MA(j))/assumptions.MA(j));
                            c_var = (assumptions.sigc^2/M)*(1+assumptions.q*(M-assumptions.MA(j))/assumptions.MA(j));
                            thetac = c_var/c_mean;
                            kc = c_mean^2/c_var;
                        else
                            c_mean = (assumptions.muc/M)*(1-assumptions.q);
                            c_var = (assumptions.sigc^2/M)*(1-assumptions.q);
                            thetac = c_var/c_mean;
                            kc = c_mean^2/c_var;
                        end
                        
                        if c_mean ~= 0
                            FiTj = gamrnd(kc,thetac,assumptions.SA(i),assumptions.MA(j));
                        else
                            FiTj = zeros(assumptions.SA(i),assumptions.MA(j));
                        end
                        Fi = [Fi,FiTj];
                    end
                    c = [c;Fi];
                end

                if assumptions.Sgen > 0
                    % Set distributions for generalists
                    % See Eq(5) in MiCRM paper
                    c_mean = assumptions.muc/M;
                    c_var = assumptions.sigc^2/M;
                    thetac = c_var/c_mean;
                    kc = c_mean^2/c_var;
                    Gen = [];

                    % Sampling for generalists per resource type
                    for j = 1:length(assumptions.MA)
                        GenTj = gamrnd(kc,thetac,assumptions.Sgen,assumptions.MA(j));
                        Gen = [Gen,GenTj];
                    end
                    c = [c;Gen];
                end
            end

            % Sampling of D

            Mw = assumptions.MA(length(assumptions.MA)); ...# waste resources
                Tw_index = M-assumptions.MA(end)+1:M; ...resources in waste
                counter = 1;

            for i = 1:length(assumptions.MA)-1 ...for all non-waste input resource

                Mi = assumptions.MA(i); ...# type i resources
                    Ti_index = counter:counter+Mi-1; ...resources in type i

                p = ones(M,1)*(1-assumptions.fs-assumptions.fw)/(M-Mi-Mw); ...background secretion
                    p(Ti_index) = assumptions.fs/Mi; ...self-secretion for type i
                    p(Tw_index) = assumptions.fw/Mw; ...waste secretion

                D(Ti_index,1:M) = drchrnd(transpose(p/assumptions.sparsity),Mi);
                counter = counter+Mi;
            end

            p = ones(M,1)*(1-assumptions.fs-assumptions.fw)/(M-Mw); ...background secretion
                p(Tw_index) = (assumptions.fs+assumptions.fw)/Mw; ...self-secretion for waste
                D(Tw_index,1:M) = drchrnd((p/assumptions.sparsity).',assumptions.MA(end));

            D = D.'; ...column is input, row is output

        end

        function params = MakeParams(assumptions)

            import essential_functions.struct2array
            import MarsBase.*

            % Generates initial resource array, consumer and metabolic matrices
            [~,R0] = MakeInitialState(assumptions);
            [c,D] = MakeMatrices(assumptions);

            % Prepare variables
            well_index = 1:assumptions.n_wells;
            params = struct();
            R0 = struct2array(R0);
            M = sum(assumptions.MA); ...total resource count
                S_tot = sum(assumptions.SA)+assumptions.Sgen; ...total species count

            % Creates struct of parameters for each well
            for i = well_index

                well_params = struct('c',c, ...sets default parameters
                    'm',ones(1,S_tot), ...zero-growth energy uptake per species
                    'w',ones(1,M), ...energy value per resource
                    'D',D, ...metabolic matrix
                    'g',ones(1,S_tot), ...growth rate/energy value per species
                    'l',zeros(1,M), ...leakage fraction per resource
                    'R0',R0(i,:), ...initial resource concentrations
                    'tau',1, ...
                    'r',ones(1,M), ...self-renewing resource supply rate
                    'sigma_max',1, ...
                    'nreg',10, ...Hill's coefficient for metabolic regulation
                    'n',2, ...Hill's coefficient for functional response
                    'extant_consumers',true(1,S_tot), ...
                    'extant_resources',true(1,M) ...
                    );

                params.('W'+string(i)) = well_params;
            end

            % To change any parameter below, set it in assumptions
            default_params = ['m' 'w' 'g' 'l' "tau" 'r' "sigma_max" 'n' "nreg"];

            for i = 1:length(default_params)
                if isfield(assumptions,default_params(i))
                    for j = well_index
                        params.('W'+string(j)).(default_params(i)) = assumptions.(default_params(i));
                    end
                end
            end
        end

        function dRdt = MakeResourceDynamics(assumptions)

            import essential_functions.*

            % Inputs:   N - 1xS_tot population array
            %           R - 1xM concentration array
            % Output:   function: N,R,params -> 1xM concentration array

            % Set up options for functional response, metabolic regulation, supply
            sigma = struct('typeI',@(R,params) (ColumnScale(params.c,R)), ...type I functional response
                'typeII',@(R,params) (ColumnScale(params.c,R)./(1+ColumnScale(params.c,R)./params.sigma_max)), ...type II ,sigma_max is 1/(handling time)
                'typeIII',@(R,params) (ColumnScale(params.c,R).^params.n./(1+(ColumnScale(params.c,R).^params.n)./params.sigma_max)) ...type III
                );

            u = struct('independent',@(x,params) (1), ...
                'energy',@(x,params) (nan20(ColumnScale((ColumnScale(x,params.w).^params.nreg).',...
                1./sum(ColumnScale(x,params.w).^params.nreg,2).').')),...
                'mass',@(x,params) (nan20(ColumnScale((x.^params.nreg).',1./sum(x.^params.nreg,2).').')) ...
                );

            h = struct('off',@(R,params) (0), ...
                'external',@(R,params) ((params.R0-R)./params.tau), ...
                'self_renewing',@(R,params) (params.r .*R .*(params.R0-R)));


            % Set up resource rate equation

            J_in = @(R,params) (ColumnScale(u.(assumptions.regulation)(ColumnScale(params.c,R),params),params.w)...
                .* sigma.(assumptions.response)(R,params));

            J_out = @(R,params) (ColumnScale(J_in(R,params),params.l) *(params.D.'));

            dRdt = @(N,R,params) ((h.(assumptions.supply)(R,params).'...
                -(ColumnScale(J_in(R,params),1./params.w).' *(N.')) + (ColumnScale(J_out(R,params),1./params.w).' *(N.'))).'...
                );
        end

        function dNdt = MakeConsumerDynamics(assumptions)

            import essential_functions.*

            % Inputs:   N - 1xS_tot population array
            %           R - 1xM concentration array
            % Output:   function: N,R,params -> 1xS_tot population array


            % Set up options for functional response, metabolic regulation, supply
            sigma = struct('typeI',@(R,params) (ColumnScale(params.c,R)), ...type I functional response
                'typeII',@(R,params) (ColumnScale(params.c,R)./(1+ColumnScale(params.c,R)./params.sigma_max)), ...type II ,sigma_max is 1/(handling time)
                'typeIII',@(R,params) (ColumnScale(params.c,R).^params.n./(1+(ColumnScale(params.c,R).^params.n)./params.sigma_max)) ...type III
                );

            u = struct('independent',@(x,params) (1), ...
                'energy',@(x,params) (nan20(ColumnScale((ColumnScale(x,params.w).^params.nreg).',...
                1./sum(ColumnScale(x,params.w).^params.nreg,2).').')),...
                'mass',@(x,params) (nan20(ColumnScale((x.^params.nreg).',1./sum(x.^params.nreg,2).').')) ...
                );

            % Set up consumer rate equation

            J_in = @(R,params) (ColumnScale(u.(assumptions.regulation)(ColumnScale(params.c,R),params),params.w)...
                .* sigma.(assumptions.response)(R,params));

            J_growth = @(R,params) (ColumnScale(J_in(R,params),(1-params.l)));

            dNdt = @(N,R,params) (params.g .*N .*((sum(J_growth(R,params),2)-params.m.').'));

        end

        function inv_info = MakeInvader(M,sampling,m,g)
            
            import MarsBase.*

            % sampling - struct:    'sampling' = sampling method (Binary, Gamma, Gaussian)
            %                       include any params needed above
            %                            binary   -> muc, q, c0, c1, 
            %                            gamma    -> muc, q, sigc
            %                            gaussian -> muc, q, sigc
            
            a_inv = struct('MA',ones(1,M), ...
                        'SA',[0],...
                        'Sgen',1,...
                        'sampling',sampling.sampling,...
                        'muc',sampling.muc,...
                        'q',sampling.q,...
                        'fs',1,...
                        'fw',1,...
                        'sparsity',1);

            if strcmp(sampling.sampling,'Binary')
                a_inv.c0 = sampling.c0;
                a_inv.c1 = sampling.c1;
            else
                a_inv.sigc = sampling.sigc;
            end

            [c,~] = MakeMatrices(a_inv);
            inv_info = struct(c=c, m=m, g=g);
        end
    
    
    end
end













