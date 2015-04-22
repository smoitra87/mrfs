function [a,b] = dummy(a,b, varargin)

fprintf('nargin : %d\n',nargin)


if length(varargin) > 0
    for i = 1:length(varargin)
        fprintf('varargin{%d} : %s\n',i,varargin{i})
    end
end

