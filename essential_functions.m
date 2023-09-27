
classdef essential_functions
    methods(Static)
        function r = drchrnd(a,n)
            % take a sample from a dirichlet distribution
            p = length(a);
            r = gamrnd(repmat(a,n,1),1,n,p);
            r = r ./ repmat(sum(r,2),1,p);
        end
        
        function a = struct2array(s)
            c = struct2cell(s);
            a = [c{:}];
            
            if istable(a)
                a = table2array(a);
            end
        end
        
        function sm = ColumnScale(m,v)
            
            % Scales each column of m by corresponding elt in v
            
            [a,b] = size(m);
            assert(b == size(v,2) | b == size(v,1),'Size of vector must equal number of columns in matrix')
            
            sm = zeros(a,b);
            
            for i = 1:b
                sm(:,i) = m(:,i).*v(i);
            end
        end
        
        function m = nan20(n)
            n(isnan(n)) = 0;
            m = n;
        end
        
        function y = ConcatNR(N,R)
            if isstruct(N) || isstruct(R)
                N = essential_functions.struct2array(N);
                R = essential_functions.struct2array(R);
            end
            y = [N,R];
        end
        
        function params_comp = CompressParams(extant_consumers,extant_resources,params,S,M)
            %{
        NOTE: extant_resources should be true for all resources if
        renewable, ie crossfeeding

        extant_consumers - 1xS_tot boolean array, true for positive
            population per species
        extand_resources - 1xM boolean array, true for positive
            concentration per resource
            %}
            
            params_comp = params;
            
            % Variables of size SxM
            assert(isequal(size(params.c),[S M]),'Size of c is not SxM');
            params_comp.c = params.c(extant_consumers,:);
            params_comp.c = params_comp.c(:,extant_resources);
           
            
            % Variables of size MxM
            assert(isequal(size(params.D),[M M]),'Size of D is not MxM');
            params_comp.D = params.D(extant_resources,:);
            params_comp.D = params_comp.D(:,extant_resources);
            
            % Variables of size 1xS
            assert(isequal(size(params.m),[1 S]),'Size of m is not 1xS');
            assert(isequal(size(params.g),[1 S]),'Size of g is not 1xS');
            params_comp.m = params.m(extant_consumers);
            params_comp.g = params.g(extant_consumers);
            
            % Variables of size 1xM
            assert(isequal(size(params.w),[1 M]),'Size of w is not 1xM');
            assert(isequal(size(params.r),[1 M]),'Size of r is not 1xM');
            assert(isequal(size(params.tau),[1 M]),'Size of tau is not 1xM');
            assert(isequal(size(params.R0),[1 M]),'Size of R0 is not 1xM');
            assert(isequal(size(params.l),[1 M]),'Size of l is not 1xM');
            params_comp.w = params.w(extant_resources);
            params_comp.r = params.r(extant_resources);
            params_comp.tau = params.tau(extant_resources);
            params_comp.R0 = params.R0(extant_resources);
            params_comp.l = params.l(extant_resources);
            
            % Consumer and resource indices
            params_comp.extant_consumers = extant_consumers;
            params_comp.extant_resources = extant_resources;
        end
        
    end
end


