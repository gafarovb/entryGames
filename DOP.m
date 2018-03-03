classdef DOP
    %DOP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        order = 2;
        coeffVector; %alpha
        x;
        bernFromSP;
        binom;
    end
    
    methods
        
        function obj = DOP(x)
            
            % obj.order = order;
            obj.bernFromSP = obj.BernsteinFromSP( obj.order);
            obj.binom = obj.nchoosekMat( obj.order);
            
            obj.x = x;
            nCoeff = (obj.order+1)^2;
            obj.coeffVector =  0.1* ones(nCoeff,1);
        end
        
        
        
        function basPoly = bernstein_i(obj,i,v)
            
            n = obj.order;
            
            basPoly  =  nchoosek(n,i)  * (v(:,1).^i).*((1-v(:,1)).^(n-i)) ;
            
        end
        
        function fun = sieve(obj,v)
            n = obj.order;
            fun = 0;
            lx =1;
            for jx = 0:n
                for ix = 0:n
                    fun =   fun  + obj.coeffVector(lx)* obj.bernstein_i(ix,v(:,1)) * obj.bernstein_i(jx,v(:,2));
                    lx = lx + 1;
                end
            end
        end
        
        
        function uniIntegral = uniIntegrals(obj,theta )
            beta1 = theta(1:2) ;
            delta1 = theta(3) ;
            beta2 = theta(4:5) ;
            delta2 = theta(6) ;
            x = obj.x;
            xB_1 = [1 x(1)]* beta1;
            xB_2 = [1 x(2)]* beta2;
            
            n = obj.order;
            
            SPfromI = obj.bernFromSP * obj.SPfromI_ix(obj.binom,delta1,xB_1);
            SPfromI(:,:,2) = obj.bernFromSP * obj.SPfromI_ix(obj.binom,delta2,xB_2);
            
            uniIntegral(:,1) = SPfromI(:,:,1)*(obj.stdIntegral(n, delta1-xB_1) - obj.stdIntegral(n, -xB_1));
            uniIntegral(:,2) = SPfromI(:,:,2)*(obj.stdIntegral(n, delta1-xB_2) - obj.stdIntegral(n, -xB_2));
            
        end
          function vecTesnorInt = vecTensorIntegrals(obj,theta )
              uniIntegr = uniIntegrals(obj,theta );
              tensorInt = kron(uniIntegr(:,1),uniIntegr(:,2)');
              
              n = size(tensorInt,1);
              vecTesnorInt = reshape(tensorInt, n*n,1);
              
          end
         
          
          function answ = evalIntegral(obj,theta)
              answ = vecTensorIntegrals(obj,theta )'*obj.coeffVector;
          end
    end
    
  
    methods (Static)
        
        function binom = nchoosekMat(order)
            
            binom  = zeros(order+1, order+1);
            for ix = 0:order
                for jx = 0:ix
                    binom(ix+1, jx+1) = nchoosek(ix,jx) ;
                end
            end
        end
        
        function conversion = BernsteinFromSP(order)
            
             
            ix =  order;
            conversion = zeros(ix+1, order+1);
            
            lineNum = 0;
            for kx = 0:ix
                lineNum = lineNum+1;
                C_i_k = nchoosek(ix,kx);
                for jx= 0:(ix-kx)
                    conversion(lineNum , 1 + kx + jx) = C_i_k * ((-1)^jx)*nchoosek(ix-kx,jx);
                end
            end
             
             
            
            
            
        end
        
        
        function  conversion = SPfromI_ix(binom,delta,xB)
            order = size(binom,1)-1;
            xbPowers = (delta-xB).^(0:order);
            
            rep=  repmat(xbPowers',1,order+1) ;
            for ix = 1:order
                rep(:,ix+1)=[0; rep(1:order,ix)];
            end
            
            deltaPowers = (delta).^(-(0:order));
            repDelta = repmat(deltaPowers',1,order+1) ;
            conversion = rep.* binom.*repDelta;

        end
        
        
        function Integr = stdIntegral(order, x)
            Integr = zeros(order+1,2);
            
            for ix =0:order
                switch ix
                    case 0
                        Integr(ix+1,:) = [0 -1];
                    case 1
                        Integr(ix+1,:) = [-1 0];
                    case 2
                        Integr(ix+1,:) = [-x 1];
                    otherwise
                        Integr(ix+1,:) = [(ix-1)*Integr(ix-1,1)- x^(ix-1)   (ix-1)*Integr(ix-1,2)];
                end
            end
            phi=[normpdf(x) ; normcdf(x) ];
            Integr = Integr * phi;
            
            
        end
    end
    
    
    
end
