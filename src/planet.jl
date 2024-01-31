module planet

include("core_module.jl")

struct Planet_object
	spice_name::Any
	spice_id::Any

	keplerian_data::Any

	name::String
end

function spice(target, observer)
	name, id
	if typeof(target) == String
		name = target
		id = SPICE.bodn2c(target)
	elseif typeof(target) == Int
		name = SPICE.bodc2s(target)
		id = target
	end
	Planet_object(name, id, nothing, name)
end

end
