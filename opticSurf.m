classdef opticSurf < handle
    properties
        n_out
        n_in
        px
        samples
        thickness
        z
        amp
        lambda
    end
   
    
    methods
        function obj = opticSurf()
        end
        
        function U_out = transmitField(obj,u)            
            U_out = 2*pi/obj.lambda * (obj.n_out-obj.n_in) * obj.z .* obj.amp .* u;  %Multiply complex field by phase and amplitude of surface
        end
        function U_prop = propagateInput(obj,u)
            U_out = obj.transmitField(u);
            [Fx,Fy] = obj.makeFreq;
            U_prop = propagate2(U_out,obj.lambda,obj.thickness,Fx,Fy);
        end
        function [Fx, Fy] = makeFreq(obj)
            fx = [1:obj.samples(2)]*obj.px;
            fy = [1:obj.samples(1)]*obj.px;
            [Fy,Fx] = meshgrid(fy,fx);
        end
    end
end