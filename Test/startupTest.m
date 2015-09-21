function rootpath = startupTest()
% initialization function for Test methods
rootpath = fileparts(pwd);
addpath(rootpath);
startup(rootpath);
end