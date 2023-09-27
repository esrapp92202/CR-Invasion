
classdef Community < handle
    
    properties
        N % array of well current species populations
        R % array of well current resource concentrations
        R0 % initial resource concentrations (1xM array)
        dNdt % consumer dynamics
        dRdt % resource dynamics
        params % struct of well params
        n_wells % number of wells
        S % number of species
        M % number of resources
        Nt % struct of well consumer trajectories
        Rt % struct of well resource trajectories
        SS % vector of species composition in initial steady state
    end
    
    methods
        
        function obj = Community(init_state,params,dynamics)
            
            % init_state - {N0,R0} where N0,R0 structs
            % params - struct of parameters for each well
            % dynamics - {dNdt,dRdt}
            
            import essential_functions.struct2array
            
            [obj.N,obj.R] = deal(struct2array(init_state{1}),struct2array(init_state{2}));
            [obj.dNdt,obj.dRdt] = deal(dynamics{1},dynamics{2});
            obj.R0 = obj.R;
            obj.params = params;
            [obj.n_wells,obj.S] = size(obj.N);
            obj.M = size(obj.R,2);

            obj.Nt = struct();
            obj.Rt = struct();

            for i = 1:obj.n_wells
                obj.Nt.('W'+string(i)) = obj.N(i,:);
                obj.Rt.('W'+string(i)) = obj.R(i,:);
            end
        end
        
        function obj = Reset(obj,init_state) % Reset state while maintaining params and dynamics

            import essential_functions.struct2array

            [obj.N,obj.R] = deal(struct2array(init_state{1}),struct2array(init_state{2}));
            obj.R0 = obj.R;
            [obj.n_wells,obj.S] = size(obj.N);
            obj.M = size(obj.R,2);

            % Remove any additional species from params
            for i = 1:obj.n_wells
                obj.params.('W'+string(i)).c = obj.params.('W'+string(i)).c(1:obj.S,:);
                obj.params.('W'+string(i)).m = obj.params.('W'+string(i)).m(1:obj.S);
                obj.params.('W'+string(i)).g = obj.params.('W'+string(i)).g(1:obj.S);
            end
            
            % Reset trajectories
            obj.Nt = struct();
            obj.Rt = struct();

            for i = 1:obj.n_wells
                obj.Nt.('W'+string(i)) = obj.N(i,:);
                obj.Rt.('W'+string(i)) = obj.R(i,:);
            end
        end
        
        function [NR_comp,params_comp,dydt,S_comp] = PrepareWell(obj,well_index,compress_consumers,compress_resources)
            
            % Compression currently nonfunctional
            
            import essential_functions.*
            
            NR = [obj.N obj.R]; % Combines N and R into single array
            NR = NR(well_index,:); % Creates NR vector for indicated well
            
            params_uncomp = obj.params.('W'+string(well_index)); % Uncompressed parameters for indicated well
            
            % Compresses parameters for large simulations
            %{
            
            % Generates boolean arrays of extant species and resources
            if compress_consumers
                extant_consumers = NR(1:obj.S)>0;
            else
                extant_consumers = true(1,obj.S);
            end
            
            if compress_resources
                extant_resources = NR(obj.S+1:end)>0;
            else
                extant_resources = true(1,size(NR,2)-obj.S)>0;
            end
            %}
            
            % While compression nonfunctional, consider all extant
            extant_consumers = true(1,obj.S);
            extant_resources = true(1,obj.M);
            
            % Compresses params, records new NR and S
            params_comp = CompressParams(extant_consumers,extant_resources,params_uncomp,obj.S,obj.M);
            S_comp = sum(extant_consumers);
            extant_NR = [extant_consumers,extant_resources];
            NR_comp = NR(extant_NR);
            
            % Generates dynamical equation using compressed params
            dN = @(NR) (obj.dNdt(NR(1:S_comp).',NR(S_comp+1:end).',params_comp));
            dR = @(NR) (obj.dRdt(NR(1:S_comp).',NR(S_comp+1:end).',params_comp));
            dydt = @(NR) ([dN(NR),dR(NR)]);
        end
        
        function [t,N_traj,R_traj] = Propagate(obj,T,dt,compress_consumers,compress_resources)
            
            t_scale = linspace(0,T,T/dt);
            N_traj = struct();
            R_traj = struct();
            
            
            for i = 1:obj.n_wells
                [NR_comp,params_comp,dydt,S_comp] = obj.PrepareWell(i,false,false);...compress_consumers,compress_resources);
                [t,NR] = ode45(@(t,NR_comp) dydt(NR_comp).',t_scale,NR_comp);
                N_traj.('W'+string(i)) = NR(:,1:S_comp);
                R_traj.('W'+string(i)) = NR(:,S_comp+1:end);
                
                obj.N(i,:) = NR(end,1:S_comp);
                obj.R(i,:) = NR(end,S_comp+1:end);
                
                obj.Nt.('W'+string(i)) = [obj.Nt.('W'+string(i));N_traj.('W'+string(i))];
                obj.Rt.('W'+string(i)) = [obj.Rt.('W'+string(i));R_traj.('W'+string(i))];
                % Compresses parameters (nonfunctional)
                %{
                
                % Updates well params to params_comp
                params_uncomp = obj.params.('W'+string(i));
                obj.params.('W'+string(i)) = params_comp;
                
                % Updates N and R
                
                
                extant_consumer_indices = find(params_comp.extant_consumers);
                extant_resource_indices = find(params_comp.extant_resources);
                
                for j = 1:sum(params_comp.extant_consumers)
                    obj.N(i,extant_consumer_indices(j)) = N_traj.('W'+string(i))(end,j);
                end
                
                for j = 1:sum(params_comp.extant_resources)
                    obj.R(i,extant_resource_indices(j)) = R_traj.('W'+string(i))(end,j);
                end
                
                
                obj.params.('W'+string(i)) = params_uncomp; % TEMPORARY
                
                %}
            end
            

        end
        
        function obj = Passage(obj,tm,ext_thres,refresh_resources)
            
            % tm - transfer matrix (A,B) -> old well A to new well B
            % refresh_resources - if true, resources are replenished by R0
            
            assert(isequal(size(tm),[obj.n_wells obj.n_wells]),'Invalid transfer matrix')
            
            % Set any negative values to 0
            obj.N = max(obj.N,0);
            obj.R = max(obj.R,0);
            
            obj.N = (obj.N.' * tm).';

            extinct_consumers = (obj.N/sum(obj.N))<ext_thres;
            obj.N(extinct_consumers) = 0;
            
            if refresh_resources
                obj.R = (obj.R.' * tm).' + obj.R0;
            else
                obj.R = (obj.R.' * tm).';
            end
        end
        
        function [Nf,Rf] = FindSteadyState(obj,thres,dt,max_iter,cons,extinct_thres)
            
            % thres - threshhold of error between propagations, (0<thres<1)
            % dt - time for each propagation step
            % max_iter - maximum iterations of propagation (int)
            % cons - # times in a row error must remain below thres (int)
            
            i = 0; % iteration counter
            j = 0; % consistency counter
            
            while j < cons && i < max_iter
                Ncomp0 = obj.N ./sum(obj.N,2);
                Rcomp0 = obj.R ./sum(obj.R,2);
                comp0 = [Ncomp0 Rcomp0];
                
                obj.Propagate(dt,1,false,false);
                
                Ncompf = obj.N ./sum(obj.N,2);
                Rcompf = obj.R ./sum(obj.R,2);
                compf = [Ncompf Rcompf];
                
                diff_array = abs(compf-comp0);
                diff = max(diff_array,[],"all");
                
                i = i + 1;
                
                if diff > thres
                    j = 0;
                else
                    j = j + 1;
                end
                
                if max(obj.N,[],"all") < extinct_thres
                    disp("Community collapsed")
                    j = cons + 1;
                end
                
            end
            
            if i == max_iter
                Nf = NaN(size(obj.N));
                Rf = NaN(size(obj.R));
            else
                Nf = obj.N;
                Rf = obj.R;
            end
            
        end

        function obj = Invade(obj,inv_info,load)
            % inv_info - struct:    c = consumption preferences of invader (1xM array)
            %                       m = zero growth rate (number)
            %                       g = basal growth rate
            % load - percent of population which invader will comprise (between 0 and 1)

            N_inv = zeros(obj.n_wells,1);

            for i = 1:obj.n_wells
                N_inv(i) = load*sum(obj.N(i,:));
                obj.params.('W'+string(i)).c = [obj.params.('W'+string(i)).c; inv_info.c];
                obj.params.('W'+string(i)).m = [obj.params.('W'+string(i)).m inv_info.m];
                obj.params.('W'+string(i)).g = [obj.params.('W'+string(i)).g inv_info.g];
                obj.params.('W'+string(i)).extant_consumers = [obj.params.('W'+string(i)).extant_consumers 1];
                obj.Nt.('W'+string(i)) = [obj.Nt.('W'+string(i)) zeros(size(obj.Nt.('W'+string(i)),1),1)];
            end

            obj.N = [obj.N N_inv];
            obj.S = obj.S + 1;
        end
    
    end
end