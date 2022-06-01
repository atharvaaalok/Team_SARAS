classdef standardAtmosphere 
% class for a standard atmosphere with only contant temperature or gradient layers
% initialized with ISA for Earth

    properties
        SL          % sea-level properties: pressure, temperature
        layer       % struct with fields startAltitude, pauseAltitude, lapseRate
                    % initialized with troposphere and stratosphere upto 32 km
        constant    % atmospheric constants: acceleration due to gravity, gas constant,
                    % ratio of specific heats, radius of Earth
    end
    
    methods
        
        % class constructor
        % initializes with ISA constants and sea level propeties
        % layers up to 32 km
        function atmos = standardAtmosphere()
            atmos.constant = standardAtmosphere.setAtmosphereConstants();
            atmos.SL = standardAtmosphere.setSeaLevelProperties();
            atmos.layer = standardAtmosphere.setLayerProperties(0,11,-6.5);
            atmos.layer(2) = standardAtmosphere.setLayerProperties(11,20,0);
            atmos.layer(3) = standardAtmosphere.setLayerProperties(20,32,1);
        end
        
        % compute temperature at a given altitude
        function T = temperature(atmos,h)
            nLayer = length(atmos.layer);
            iLayer = 0;
            isComputed = 0;
            while (iLayer < nLayer) && ~isComputed
                iLayer = iLayer + 1;
                if h <= atmos.layer(iLayer).pauseAltitude
                    if iLayer == 1
                        h0 = 0;
                        T0 = atmos.SL.temperature;
                    else
                        h0 = atmos.layer(iLayer).startAltitude;
                        T0 = atmos.temperature(h0);
                    end
                    lambda = atmos.layer(iLayer).lapseRate;
                    T = T0 + (h-h0)*lambda;
                    isComputed = 1;
                end
            end
        end
        
        % compute pressure at a given altitude
        function p = pressure(atmos,h)
            gSL = atmos.constant.gravity;
            R = atmos.constant.gasConstant;
            nLayer = length(atmos.layer);
            iLayer = 0;
            isComputed = 0;
            while (iLayer < nLayer) && ~isComputed
                iLayer = iLayer + 1;
                if h <= atmos.layer(iLayer).pauseAltitude
                    if iLayer == 1
                        h0 = 0;
                        p0 = atmos.SL.pressure;
                        T0 = atmos.SL.temperature;
                    else
                        h0 = atmos.layer(iLayer).startAltitude;
                        p0 = atmos.pressure(h0);
                        T0 = atmos.temperature(h0);
                    end
                    lambda = atmos.layer(iLayer).lapseRate;
                    if lambda
                        p = p0*(atmos.temperature(h)./T0).^(-gSL./lambda./R*1000);
                    else
                        p = exp(-gSL/R/T0*(h-h0)*1000)*p0;
                    end
                    isComputed = 1;
                end
            end
        end
        
        % compute density at a given altitude
        function rho = density(atmos,h)
            R = atmos.constant.gasConstant;
            rho = atmos.pressure(h)./R./atmos.temperature(h);
        end
        
        % compute viscosity at a given altitude
        function visc = viscosity(atmos,h)
            T = atmos.temperature(h);
            visc = atmos.computeViscocity(T);
        end
        
        % compute speed of sound given the ambient temperature
        function a = computeSoundSpeed(atmos,T)
            gamma = atmos.constant.specificHeatRatio;
            R = atmos.constant.gasConstant;
            a = sqrt(gamma*R*T);
        end
        
        % compute speed of sound at a given altitude
        function a = soundSpeed(atmos,h)
            T = atmos.temperature(h);
            a = atmos.computeSoundSpeed(T);
        end
        
        % add a new layer specifying start and pause altitudes
        % and lapse rate in it
        function atmos = addLayer(atmos,startAltitude, pauseAltitude, lapseRate)
            iLayer = length(atmos.layer) + 1;
            atmos(iLayer).startAltitude = startAltitude;
            atmos(iLayer).pauseAltitude = pauseAltitude;
            atmos(iLayer).lapseRate = lapseRate;
        end
        
        % convert geometric to geopotential altitude
        function h = geopotentialAltitude(atmos,hGeom)
            RE = atmos.constant.radius;
            h = hGeom.*(RE./(RE+hGeom));
        end
             
    end  
        
   
    methods(Static)
        
        % structure with sea level pressure and temperature corresponding to ISA
        % used during instance creation
        function propsSL = setSeaLevelProperties()
            propsSL.pressure = 101325;
            propsSL.temperature = 288.15;
        end
        
        % Create a structure with ISA constants
        % used during object construction
        function atmosConstant = setAtmosphereConstants()
            atmosConstant.gravity = 9.80665;
            atmosConstant.gasConstant = 287.05287;
            atmosConstant.specificHeatRatio = 1.4;
            atmosConstant.radius = 6356.766;
        end  
        
        % Create a structure with layer properties
        % used in the constructor method of the class
        function layerProps = setLayerProperties(startAltitude, pauseAltitude, lapseRate)
            layerProps.startAltitude = startAltitude;
            layerProps.pauseAltitude = pauseAltitude;
            layerProps.lapseRate = lapseRate;
        end
        
        % Compute air viscosity for a given temperature
        % Sutherland's law
        function viscosity = computeViscocity(T)
            viscosity = 1.458e-6*(T).^(1.5)./(T+110.4);
        end
        
        % Convert degrees to Kelvin
        function T = degree2Kelvin(d)
            T = 273.15 + d;
        end
        
    end
end

