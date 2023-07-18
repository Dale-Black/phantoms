### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 535aeec6-2037-11ee-2453-e5c6060b1f95
# ╠═╡ show_logs = false
begin
	using Pkg; Pkg.activate(temp = true)
	Pkg.add(["ImagePhantoms", "MIRTjim", "Unitful", "ImageGeoms", "CairoMakie", "PlutoUI", "Sinograms", "AxisArrays", "DataFrames", "CSV", "Interpolations"])
	Pkg.add(url = "https://github.com/Dale-Black/Attenuations.jl")

	using ImagePhantoms
	using MIRTjim
	using Unitful: mm, cm, unit, °, ustrip, g
	using ImageGeoms: ImageGeom, MaskCircle
	using CairoMakie
	using PlutoUI
	using Sinograms
	using AxisArrays
	using DataFrames
	using CSV
	using Interpolations
	using Attenuations
end

# ╔═╡ 2282b1df-c798-47c9-8e8f-01282dac6de0
TableOfContents()

# ╔═╡ 8b0a6f89-4305-4e55-8f22-49d23629be95
md"""
# Prepare Image & CT Geometry
"""

# ╔═╡ 75b81c11-96a7-41de-a911-f87b0466ef20
md"""
## Image Geometry
"""

# ╔═╡ 67edb5cf-dda7-4cac-888c-f0c1648a8612
dims = (300, 300, 20)

# ╔═╡ 3a322442-93bb-447d-a672-831a075d4272
deltas = (1mm, 1mm, 1mm);

# ╔═╡ de15c717-7d9b-4b47-801a-7e8cf0352cdb
im_geom = ImageGeom(;dims = dims, deltas = deltas);

# ╔═╡ 1bb31f52-1037-43d7-91a3-6f4f818971f8
md"""
## CT Geometry
"""

# ╔═╡ a8ab1a0b-7729-43b9-9044-d5238846cdf5
md"""
**Sinograms.jl**

CtFanArc() Fields [helpful link](https://github.com/ismrmrd/ismrmrd-paper/blob/master/code/extern/irt/fbp/arch/ct_geom.m)
- Detector
  - `ns`: number of horizontal samples
  - `nt`: number of vertical samples
  - `ds`: horizontal sample spacing
  - `dt`: vertical sample spacing
  - `offset_s`: unitless detector center offset (usually 0 or 0.25)
  - `offset_t`: unitless, usually 0
- Source
  - `na`: number of angular samples
  - `orbit`: source orbit in degrees (or Unitful), 360 for fan beam, 180 for parallel
  - `orbit_start`: starting angle in degrees (or Unitful) orbit and orbit_start must both be unitless (degrees) or have same units.
  - `source_offset`: usually 0
  - `dsd`: distance from source to detector
  - `dod`: distance from origin to detector
  - `src::CtSource`: describes the X-ray CT source trajectory. Primary support for CtSourceCircle().

Units:
- `ds`, `dt`, `source_offset`, `dsd`, `dod` must all be unitless or have the same units.

**[Canon Aquilion One Notes](https://link.springer.com/article/10.1007/s10554-020-01961-y/tables/1)**
- X-ray source filtration: 7.5 mm Al
- Detector pixel pitch: 0.5 mm
- Number of projection angles: 900
- Source to object distance: 60.0 cm
- Object to detector distance: 47.2 cm
- Focal spot size: 0.9 mm
"""

# ╔═╡ 14891577-24a3-4507-9e42-1c0c723e11f0
# ps_canon_aquilion_one_estimate = (
# 	ns = 888*2,
# 	nt = 64,
# 	ds = 0.5mm,
# 	dt = 0.5mm,
# 	offset_s = 1.25,
# 	offset_t = 0.0,
# 	na = 984,
# 	orbit = 360.0,
# 	orbit_start = 0.0,
# 	source_offset = 0.0mm,
# 	dod = 472mm,
# 	dsd = 600mm + 472mm
# )

# ╔═╡ 432d5f93-70c1-4985-9f19-b066f4dc77a0
md"""
**GE Lightspeed System Geometry**
doi: 10.1109/TMI.2006.882141
"""

# ╔═╡ 61e34bbc-232f-499d-9343-ea944e76a760
ps_ge_lightspeed = (
	 ns = 888,
	 nt = 64,
	 ds = 1.0239mm,
	 dt = 1.0964mm,
	 offset_s = 1.25,
	 offset_t = 0.0,
	 na = 984,
	 orbit = 360.0,
	 orbit_start = 0.0,
	 source_offset = 0.0mm,
	 dsd = 949.075mm,
	 dod = 408.075mm,
)

# ╔═╡ e3ca9f1b-1543-4b92-804b-b684170cdb70
ct_geom_plot3(CtFanArc(;ps_ge_lightspeed...), im_geom)

# ╔═╡ 384b4521-9b5d-4214-ac10-77f9c2fc292f
ct_geom = CtFanArc(;ps_ge_lightspeed...);

# ╔═╡ 5a92e629-b034-4615-91c2-bc3816790898
md"""
# Prepare Phantom
"""

# ╔═╡ 4a764662-5bef-41be-9caa-15212130d4e5
md"""
## Polyenergetic Materials

The linear attenuation coefficient (often symbolized by the Greek letter μ) of a material describes how much a beam of x-ray or gamma radiation is attenuated (i.e., reduced in intensity) by that material. It is a measure of the probability per unit path length that a photon will interact with the material and lose energy.

The unit of the linear attenuation coefficient is typically inverse length (1/length), such as inverse centimeters (cm⁻¹). A higher value of the linear attenuation coefficient means the material is more effective at attenuating the radiation.

In the context of medical imaging, for example, materials like bone that have high linear attenuation coefficients appear white on an x-ray image because they absorb more of the x-ray photons, while materials like air and soft tissues that have lower coefficients appear darker because more of the x-rays pass through them.

The linear attenuation coefficient depends on the energy of the radiation, as well as the atomic number and density of the material.
"""

# ╔═╡ 1d5102a1-6e7a-44a1-ab8d-c71b40450518
md"""
### Load Spectra
"""

# ╔═╡ 4dc4addd-c371-4cdc-b688-78ec19f1fc8e
dir_spectra = "/Users/daleblack/Library/CloudStorage/GoogleDrive-djblack@uci.edu/My Drive/dev/Molloilab/phantoms/spectra"

# ╔═╡ a57e0b1c-3d45-4f97-b3ac-61aa05488b03
function xray_read_spectra_char(dir_spectra, kvp)
	MM = length(kvp)
	energies = []
	spectra = []
	for i in 1:MM
		tmp = CSV.read(joinpath(dir_spectra, "spectra_$(kvp[i]).csv"), DataFrame)
		energy = tmp[:, 1]
		# The Wilderman/Sukovic spectra must be scaled by energy!
		spectrum = tmp[:, 2] .* energy
		push!(energies, energy)
		push!(spectra, spectrum)
	end
	
	# interpolate onto same energy sampling
	MM = length(spectra);
	tmp = zeros(MM, 1)
	for i in 1:MM
		tmp[i] = maximum(energies[i])
	end
	
	en = energies[findmax(tmp)[2]]
	sp = zeros(length(en), MM)
	for mm in 1:MM
		extrap = linear_interpolation(energies[mm], spectra[mm], extrapolation_bc = Line())
		sp[:, mm] = extrap(en) # spectrums
	end
	return en, sp 
end

# ╔═╡ c9678a44-7c39-4544-aeed-c46af49a67fd
function xray_read_spectra(dir_spectra; kvp = [80, 120, 140])
	en, sp = xray_read_spectra_char(dir_spectra, kvp)
	MM = size(sp, 2) # Number of spectra
	Ide = sp .* vcat(zeros(1, MM), repeat(diff(en), 1, MM)) # (N, M) # Differential spectrum
	I = sum(Ide) # [1 M] spectrum integral
	eff_mean = en' * Ide ./ I # mean energy of the incident spectrum.
	sp_at_eff_mean = zeros(MM)
	for mm in 1:MM
		extrap = linear_interpolation(en, sp[:, mm]; extrapolation_bc = 0);
		sp_at_eff_mean[mm] = extrap(eff_mean[mm])
	end
	return en, sp, Ide, I, kvp, sp_at_eff_mean
end

# ╔═╡ b4a092e9-b04b-489c-aa6c-9e94bfac67c7
en, sp, Ide, I, kvp, sp_at_eff_mean = xray_read_spectra(dir_spectra; kvp = [80, 100, 140])

# ╔═╡ a4dd04aa-1a51-429a-908b-5bc712def263
let
	f = Figure()
	ax = CairoMakie.Axis(
		f[1, 1],
		title = "X-ray Spectra",
		xlabel = "Energy (kVp)",
		ylabel = "Photons / mAs-sr per energy bin"
	)
	for i in Base.axes(sp, 2)
		lines!(en, sp[:, i], label = "$(kvp[i]) kVp Spectrum")
	end
	axislegend(ax)
	f
end

# ╔═╡ c5236013-01a9-40cc-a9cc-0641c5028c10
md"""
### Line Integrals
"""

# ╔═╡ a07655e8-85f4-4134-b773-e8966c064eed
function de_ftab_sls(; sl = [], s_n = [], s_min = [], s_max = [])
	# number of samples of the material integrals "s"
	if isempty(s_n)
		if isempty(sl)
			s_n = [45 43]
		else
			for ll in 1:length(sl)
				s_n[1, ll] = length(sl[ll]);
			end
		end
	end

	# minimum material "integrals"
	if isempty(s_min)
		if isempty(sl)
			s_min = [0 0]
		else
			for ll in 1:length(sl)
				s_min[1, ll] = minimum(sl[ll])
			end
		end
	end

	# maximum material "integrals"
	if isempty(s_max)
		if isempty(sl)
			# soft max: 50cm * 1g/cc
			# bone max: 15cm * 2g/cc (for now)
			s_max = [50 30];
		else
			for ll in 1:length(sl)
				s_max[1, ll] = maximum(sl[ll]);
			end
		end
	end

	if isempty(sl)
	    for ll in 1:length(s_n)
	        push!(sl, range(s_min[ll], s_max[ll], length = s_n[ll]))
	    end
	end

	return sl, s_n, s_min, s_max
end

# ╔═╡ 4dff0b84-3036-4879-a04f-854c89880993
md"""
### Linear Attenuation Coefficients
"""

# ╔═╡ 2790df00-ac91-4c25-883f-5343505f8009


# ╔═╡ a27c7a79-eb7e-4018-bcf4-0fb7ac9f311e
md"""
### Convert to HUs
"""

# ╔═╡ 3e2112e3-5118-4870-9257-a8324dd924d8
μ_to_hu(μ0) = 1000 * (μ0 - μ_water) / μ_water

# ╔═╡ 7c2b12f7-08b3-4b93-a4b4-635963f793c2
md"""
## QRM Geometry & Materials
"""

# ╔═╡ 58b746d1-abce-48c7-b991-3de8db3dafc6
md"""
### Thorax
"""

# ╔═╡ b677c825-fb2c-4266-a5ac-9849bd3b9cdf
begin
	# Thorax
	center = (0mm, 0mm, 0mm)
	width = (150mm, 100mm, 15mm) # x radius, y radius, height
	angles = (0, 0, 0)
	
	thorax = cylinder(center, width, angles, soft_tissue_hus[ax])
end;

# ╔═╡ 49388590-772e-4893-a508-e366ea84664f
md"""
### Lungs
"""

# ╔═╡ f0976876-6076-47ab-8fa4-49373f42455a
begin
	# Lungs
	center_left_lung = (-82.5mm, -15mm, center[3])
	center_right_lung = (82.5mm, -15mm, center[3])
	
	width_lung = (27.5mm, 60mm, 15mm)

	angles_right_lung = (π/7, 0, 0)
	angles_left_lung = (-π/7, 0, 0)

	left_lung = cylinder(center_left_lung, width_lung, angles_right_lung, lung_hus[ax])
	right_lung = cylinder(center_right_lung, width_lung, angles_left_lung, lung_hus[ax])
end;

# ╔═╡ 5fc95c85-ffba-4fa0-bdb7-f30972a60049
md"""
### Heart
"""

# ╔═╡ e81b69cc-e86f-408b-93e0-a77e2cefab1a
begin
	# Heart
	width_heart = (40mm, 40mm, 15mm)
	
	heart = cylinder(center, width_heart, angles, myocardium_hus[ax])
end;

# ╔═╡ e2819635-c8f8-48ad-a599-f6426d90dffb
md"""
### Spine
"""

# ╔═╡ b813c1c2-c717-44a4-8b51-d9888f5b8b32
begin
	# Spine insert
	center_spine1 = (0mm, -57mm, 0mm)
	width_spine1 = (15mm, 15mm, 15mm)
	center_spine2 = (0mm, -87mm, 0mm)
	width_spine2 = (22mm, 8mm, 15mm)
	center_spine3 = (0mm, -74mm, 0mm)
	width_spine3 = (4mm, 35mm, 15mm)
	
	spine1 = cylinder(center_spine1, width_spine1, angles, bone_hus[ax])
	spine2 = cuboid(center_spine2, width_spine2, (3π/2,0,0), bone_hus[ax])
	spine3 = cuboid(center_spine3, width_spine3, (3π/2,0,0), bone_hus[ax])
end;

# ╔═╡ 5ed4bd19-9033-4a0c-82d6-7791d28f8b3b
md"""
### Coronary Artery Calcium Inserts
"""

# ╔═╡ 8dc93793-d453-4dc1-87d2-b49f8a262048
begin
	# Small Calcium Inserts (1 mm)
	width_small_insert = (0.5mm, 0.5mm, 1mm)
	
	small_insert_low_density = cylinder((5mm, 5mm, 0mm), width_small_insert, angles, insert_200_hus[ax])
	small_insert_medium_density = cylinder((-5mm, 5mm, 0mm), width_small_insert, angles, insert_400_hus[ax])
	small_insert_high_density = cylinder((0mm, -5mm, 0mm), width_small_insert, angles, insert_800_hus[ax])
end;

# ╔═╡ 48b4fd90-ed1c-43c5-9a2d-fb9ca33183f7
begin
	# Medium Calcium Inserts (1 mm)
	width_medium_insert = (1.5mm, 1.5mm, 3mm)
	
	medium_insert_low_density = cylinder((10mm, 10mm, 0mm), width_medium_insert, angles, insert_200_hus[ax])
	medium_insert_medium_density = cylinder((-10mm, 10mm, 0mm), width_medium_insert, angles, insert_400_hus[ax])
	medium_insert_high_density = cylinder((0mm, -10mm, 0mm), width_medium_insert, angles, insert_800_hus[ax])
end;

# ╔═╡ 1b374271-6806-4680-a5d0-0955945af5a6
begin
	# Large Calcium Inserts (1 mm)
	width_large_insert = (2.5mm, 2.5mm, 5mm)
	
	large_insert_low_density = cylinder((15mm, 15mm, 0mm), width_large_insert, angles, insert_200_hus[ax])
	large_insert_medium_density = cylinder((-15mm, 15mm, 0mm), width_large_insert, angles, insert_400_hus[ax])
	large_insert_high_density = cylinder((0mm, -15mm, 0mm), width_large_insert, angles, insert_800_hus[ax])
end;

# ╔═╡ 01d958a1-228f-4344-b40b-492506e210d1
objects = [
	thorax
	left_lung
	right_lung
	heart
	spine1
	spine2
	spine3
	small_insert_low_density
	small_insert_medium_density
	small_insert_high_density
	medium_insert_low_density
	medium_insert_medium_density
	medium_insert_high_density
	large_insert_low_density
	large_insert_medium_density
	large_insert_high_density
];

# ╔═╡ d152f444-44b3-4f03-baed-8ea67fbd77a7
Sinograms.axes(ig1)

# ╔═╡ 61c53fb9-867c-48b2-91f6-b32ed721d222
groundtruth_phantom = phantom(Sinograms.axes(ig1)..., objects);

# ╔═╡ 49a83445-c2f2-4d43-8dcd-343cc73fd88a
@bind z1 PlutoUI.Slider(Base.axes(groundtruth_phantom, 3); show_value = true, default = div(size(groundtruth_phantom, 3), 2))

# ╔═╡ 55d7b8e6-6330-4039-a617-238a835b0daa
let
	f = Figure(resolution=(1200, 1200))

	ax = CairoMakie.Axis(
		f[1, 1],
		title="phantom"
	)
	hm = CairoMakie.heatmap!(ustrip.(groundtruth_phantom[:, :, z1]); colormap=:grays)
	Colorbar(f[:, end+1], hm)
	hidedecorations!(ax, ticks = false, ticklabels = false)
	f
end

# ╔═╡ 387b30a6-67f6-4fe1-9fca-775473e81a30
md"""
# Simulate Scans
"""

# ╔═╡ 9c6343a7-2d5e-4548-ae16-8fd87af07d03
md"""
## Get Projection Data
Using `rays` from Sinograms.jl and `radon` from ImagePhantoms.jl
"""

# ╔═╡ 1eda2628-557a-4b26-ac5a-a5e1448f7518
# projection_data = radon(rays(scanner), objects);

# ╔═╡ fb638c1d-e734-4b21-b32b-0bc7880355e2
md"""
## Reconstruct To Image Data
Convert from projection data to image data, using filtered-back projection
"""

# ╔═╡ 6d7ca293-2f35-418e-9211-069c9e686cb3
begin
	plan = plan_fbp(scanner, ig1; window = Window(Hamming(), 1.0))
	image_data = fdk(plan, projection_data)
end;

# ╔═╡ dc86e054-9f05-4db6-aee3-69e379fdc26c
@bind z2 PlutoUI.Slider(Base.axes(image_data, 3), show_value=true, default = div(size(image_data, 3), 2))

# ╔═╡ 5c949f69-f79d-4cf0-9dd0-698ad5dd94b2
let
	f = Figure(resolution = (1400, 2000))

	ax = CairoMakie.Axis(
		f[1, 1],
		title="Ground Truth Phantom"
	)
	CairoMakie.heatmap!(ustrip.(groundtruth_phantom[:, :, z2]); colormap=:grays)
	
	ax = CairoMakie.Axis(
		f[1, 2],
		title="Reconstructed Image Data"
	)
	CairoMakie.heatmap!(ustrip.(image_data[:, :, z2]); colormap=:grays)

	ax = CairoMakie.Axis(
		f[2, 1],
		title="Projection Data"
	)
	CairoMakie.heatmap!(ustrip.(projection_data[:, :, z2]); colormap=:grays)

	ax = CairoMakie.Axis(
		f[2, 2],
		title="Error"
	)
	err = ustrip.(image_data[:, :, z2]) - ustrip.(groundtruth_phantom[:, :, z2])
	CairoMakie.heatmap!(err; colormap=:grays)
	
	f
end

# ╔═╡ Cell order:
# ╠═535aeec6-2037-11ee-2453-e5c6060b1f95
# ╠═2282b1df-c798-47c9-8e8f-01282dac6de0
# ╟─8b0a6f89-4305-4e55-8f22-49d23629be95
# ╟─75b81c11-96a7-41de-a911-f87b0466ef20
# ╠═67edb5cf-dda7-4cac-888c-f0c1648a8612
# ╠═3a322442-93bb-447d-a672-831a075d4272
# ╠═de15c717-7d9b-4b47-801a-7e8cf0352cdb
# ╟─1bb31f52-1037-43d7-91a3-6f4f818971f8
# ╟─a8ab1a0b-7729-43b9-9044-d5238846cdf5
# ╠═14891577-24a3-4507-9e42-1c0c723e11f0
# ╟─432d5f93-70c1-4985-9f19-b066f4dc77a0
# ╠═61e34bbc-232f-499d-9343-ea944e76a760
# ╠═e3ca9f1b-1543-4b92-804b-b684170cdb70
# ╠═384b4521-9b5d-4214-ac10-77f9c2fc292f
# ╟─5a92e629-b034-4615-91c2-bc3816790898
# ╟─4a764662-5bef-41be-9caa-15212130d4e5
# ╟─1d5102a1-6e7a-44a1-ab8d-c71b40450518
# ╠═4dc4addd-c371-4cdc-b688-78ec19f1fc8e
# ╠═c9678a44-7c39-4544-aeed-c46af49a67fd
# ╠═a57e0b1c-3d45-4f97-b3ac-61aa05488b03
# ╠═b4a092e9-b04b-489c-aa6c-9e94bfac67c7
# ╟─a4dd04aa-1a51-429a-908b-5bc712def263
# ╟─c5236013-01a9-40cc-a9cc-0641c5028c10
# ╠═a07655e8-85f4-4134-b773-e8966c064eed
# ╟─4dff0b84-3036-4879-a04f-854c89880993
# ╠═2790df00-ac91-4c25-883f-5343505f8009
# ╟─a27c7a79-eb7e-4018-bcf4-0fb7ac9f311e
# ╠═3e2112e3-5118-4870-9257-a8324dd924d8
# ╟─7c2b12f7-08b3-4b93-a4b4-635963f793c2
# ╟─58b746d1-abce-48c7-b991-3de8db3dafc6
# ╠═b677c825-fb2c-4266-a5ac-9849bd3b9cdf
# ╟─49388590-772e-4893-a508-e366ea84664f
# ╠═f0976876-6076-47ab-8fa4-49373f42455a
# ╟─5fc95c85-ffba-4fa0-bdb7-f30972a60049
# ╠═e81b69cc-e86f-408b-93e0-a77e2cefab1a
# ╟─e2819635-c8f8-48ad-a599-f6426d90dffb
# ╠═b813c1c2-c717-44a4-8b51-d9888f5b8b32
# ╟─5ed4bd19-9033-4a0c-82d6-7791d28f8b3b
# ╠═8dc93793-d453-4dc1-87d2-b49f8a262048
# ╠═48b4fd90-ed1c-43c5-9a2d-fb9ca33183f7
# ╠═1b374271-6806-4680-a5d0-0955945af5a6
# ╠═01d958a1-228f-4344-b40b-492506e210d1
# ╠═d152f444-44b3-4f03-baed-8ea67fbd77a7
# ╠═61c53fb9-867c-48b2-91f6-b32ed721d222
# ╟─49a83445-c2f2-4d43-8dcd-343cc73fd88a
# ╟─55d7b8e6-6330-4039-a617-238a835b0daa
# ╟─387b30a6-67f6-4fe1-9fca-775473e81a30
# ╟─9c6343a7-2d5e-4548-ae16-8fd87af07d03
# ╠═1eda2628-557a-4b26-ac5a-a5e1448f7518
# ╟─fb638c1d-e734-4b21-b32b-0bc7880355e2
# ╠═6d7ca293-2f35-418e-9211-069c9e686cb3
# ╟─dc86e054-9f05-4db6-aee3-69e379fdc26c
# ╟─5c949f69-f79d-4cf0-9dd0-698ad5dd94b2
