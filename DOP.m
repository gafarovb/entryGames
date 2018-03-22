classdef DOP<handle
    %DOP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        order = 2;
        outcome = [1;0];
        x;
        bernFromSP;
        binom;
    end
    
    methods
        
        function obj = DOP(x,outcome )
            
            % obj.order = order;
            obj.outcome = outcome;
            obj.bernFromSP = obj.BernsteinFromSP( obj.order);
            obj.binom = obj.nchoosekMat( obj.order);
            
            obj.x = x;
          end
        function basPoly = bernstein_i(obj,i,v)
            basPoly  =  nchoosek(obj.order,i)  * (v(:,1).^i).*((1-v(:,1)).^(obj.order-i)) ;
        end
        function fun = sieve(obj,coeffVector,v)
            n = obj.order;
            fun = 0;
            lx =1;
            for jx = 0:n
                for ix = 0:n
                    fun =   fun  + coeffVector(lx)* obj.bernstein_i(ix,v(:,1)) * obj.bernstein_i(jx,v(:,2));
                    lx = lx + 1;
                end
            end
        end
        function [mu,Sigma,xB,delta]= muSigma(obj,theta)
            beta1 = theta(1:2) ;
            delta1 = theta(3) ;
            beta2 = theta(4:5) ;
            delta2 = theta(6) ;
            x = obj.x;
            xB_1 = [1 x(1)]* beta1;
            xB_2 = [1 x(2)]* beta2;
            
            
            xB = [xB_1;xB_2];
            delta = [delta1;delta2];
            mu = ones(2,1) - xB./delta ;
            Sigma = [1/delta1^2 0;...
                0 1/delta2^2];
        end
        function uniIntegral = uniIntegrals(obj,theta )
            [mu,Sigma,xB,delta] = muSigma(obj,theta);
            for i = 1:2
                SPfromI(:,:,i) = obj.bernFromSP * obj.SPfromI_ix(delta(i),xB(i));
                ulim = (1 - mu(i))/sqrt(Sigma(i,i));
                llim = ulim - 1 /sqrt(Sigma(i,i));
                uniIntegral(:,i) = SPfromI(:,:,i)*(obj.stdIntegral(ulim ) - obj.stdIntegral(llim ));
            end
            
        end
        function vecTesnorInt = vecTensorIntegrals(obj,theta )
            uniIntegr = uniIntegrals(obj,theta );
            tensorInt = (uniIntegr(:,1)*uniIntegr(:,2)');
            
            n = size(tensorInt,1);
            vecTesnorInt = reshape(tensorInt, n*n,1);
            
        end
        function  conversion = SPfromI_ix(obj,delta,xB)
            binom = obj.binom;
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
        function answ = evalIntegral(obj,theta,coeffVector)
            answ = vecTensorIntegrals(obj,theta )'*coeffVector;
        end
        function Integr = stdIntegral(obj, x)
            % this function was unit tested for order up to 4 on 3/20/2018
            order = obj.order;
            Integr = zeros(order+1,2);
            
            for ix =0:order
                switch ix
                    case 0
                        Integr(ix+1,:) = [0 1];
                    case 1
                        Integr(ix+1,:) = [1 0];
                    case 2
                        Integr(ix+1,:) = [x 1];
                    otherwise
                        Integr(ix+1,:) = [(ix-1)*Integr(ix-1,1)+ x^(ix-1)   (ix-1)*Integr(ix-1,2)];
                end
            end
            phi=[-normpdf(x) ; normcdf(x) ];
            Integr = Integr * phi;
            
            
        end
        
        
        function intAy = integralAy(obj,theta)
            [mu,Sigma,xB,delta]= muSigma(obj,theta);
            
            for i = 1:2
                ulim(i) = (1 - mu(i))/sqrt(Sigma(i,i));
                llim(i) = ulim(i) - 1 /sqrt(Sigma(i,i));           
            end
            
            switch 10*obj.outcome(1) + obj.outcome(2)
                case 00 
                    intAy = prod(1-normcdf(ulim));
                case 11 
                     intAy = prod(normcdf(llim));
                case 10
                     intAy = (1-normcdf(llim(2)))*normcdf(llim(1))  + (1-normcdf(ulim(2)))*(normcdf(ulim(1)) - normcdf(llim(1)) );
                case 01
                     intAy =  (1-normcdf(llim(1)))*normcdf(llim(2))  + (1-normcdf(ulim(1)))*(normcdf(ulim(2)) - normcdf(llim(2)) );
                otherwise
                         error('Wrong value for the outcome')

            end
 
        end
        
        
        function ll = loglikelihood(obj,gamma)
            if size(obj,1)>1
               ll = 0;
               for ix = 1:size(obj,1)
                  ll = ll+ loglikelihood(obj(ix),gamma);
               end
            else
            theta = gamma(1:6);
            hMat  = reshape(gamma(7:end),(obj.order+1)^2,3);
            
            switch 10*obj.outcome(1) + obj.outcome(2)
                case 00
                    coeffVector = hMat(:,1);
                case  11
                    coeffVector = hMat(:,2);
                case 10
                    coeffVector = hMat(:,3);
                case 01
                    coeffVector = 1- sum( hMat,2);
            end
             
             
             ll =log (logBarrier( evalIntegral(obj,theta,coeffVector)+integralAy(obj,theta)));
            end
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
        
    end
    
    
    
end
