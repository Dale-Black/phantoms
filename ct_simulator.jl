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
	Pkg.add(["ImagePhantoms", "MIRTjim", "Unitful", "ImageGeoms", "CairoMakie", "PlutoUI", "Sinograms", "AxisArrays"])
	Pkg.add(url = "https://github.com/Dale-Black/Attenuations.jl")

	using ImagePhantoms
	using MIRTjim: jim
	using Unitful: mm, cm, unit, °, ustrip, g
	using ImageGeoms: ImageGeom, MaskCircle
	using CairoMakie
	using PlutoUI
	using Sinograms
	using AxisArrays
	# using FFAST
	using Attenuations
end

# ╔═╡ 2282b1df-c798-47c9-8e8f-01282dac6de0
TableOfContents()

# ╔═╡ 1cd5f63c-36af-4c72-bccd-5c8f25dff3fe
md"""
# Prepare Materials

The linear attenuation coefficient (often symbolized by the Greek letter μ) of a material describes how much a beam of x-ray or gamma radiation is attenuated (i.e., reduced in intensity) by that material. It is a measure of the probability per unit path length that a photon will interact with the material and lose energy.

The unit of the linear attenuation coefficient is typically inverse length (1/length), such as inverse centimeters (cm⁻¹). A higher value of the linear attenuation coefficient means the material is more effective at attenuating the radiation.

In the context of medical imaging, for example, materials like bone that have high linear attenuation coefficients appear white on an x-ray image because they absorb more of the x-ray photons, while materials like air and soft tissues that have lower coefficients appear darker because more of the x-rays pass through them.

The linear attenuation coefficient depends on the energy of the radiation, as well as the atomic number and density of the material.
"""

# ╔═╡ ae2e029b-0590-4ed4-911d-9d4a80dea2af
energies = [80keV, 100keV, 120keV, 135keV]

# ╔═╡ fa5406de-baa2-4f66-b883-2c551ec813b8
axs = AxisArrays.Axis{:energy}(energies)

# ╔═╡ bb1f98c9-7163-4eaa-a1a5-8e6009a53870
ax = AxisArrays.Axis{:energy}(100keV)

# ╔═╡ 1e62720c-e9bb-4eec-a7b9-e39e6731c574
md"""
## Water
"""

# ╔═╡ 1fc443ff-39a6-4f82-a582-6aa35ddd064f
water_lacs = μ(Materials.water, energies)

# ╔═╡ 6deb0f71-a009-4ac0-b48a-ec334280eaec
function convert_hu(lacs, energies)
    values = [1000 .* ((lacs[AxisArrays.Axis{:energy}(e)] ./ water_lacs[AxisArrays.Axis{:energy}(e)]) .- 1) for e in energies]
    AxisArray(values, axs)
end

# ╔═╡ 2d3c740d-da12-4f61-899e-ee5252e33b48
water_hus = convert_hu(water_lacs, energies)

# ╔═╡ 8728f157-87c1-4a5e-a376-9f71f9e3d82a
water_hus[ax]

# ╔═╡ 1f5b5c39-a125-45bd-b046-6f11060c2d78
md"""
## Air
"""

# ╔═╡ 6e3577ae-3ea6-4ce0-9f2d-381111aa05d2
air_lacs = μ(Materials.air, energies)

# ╔═╡ 9e7cae05-ac89-4ba7-a5c9-2a8e102ecb39
air_hus = convert_hu(air_lacs, energies)

# ╔═╡ 482c2bf8-1234-4e08-ac2f-7f40dda8a1f4
air_hus[ax]

# ╔═╡ 5a75e9ae-bd8e-4a15-9690-c00d95749f41
md"""
## Lung
"""

# ╔═╡ 1ed21aef-14bc-49a9-a397-c6c5d15e7c31
lung_tissue_lacs = μ(Materials.lung, energies)

# ╔═╡ 8c3dcf4f-c3d5-444e-8693-e487e1ed54e5
lung_lacs = AxisArray((0.75*air_lacs) + (0.25*lung_tissue_lacs), axs)

# ╔═╡ cc2b185d-e8b4-4283-a9d2-5af5755606d0
lung_hus = convert_hu(lung_lacs, energies)

# ╔═╡ 9a54542a-d333-4b08-ad8e-17be90b8b0fc
lung_hus[ax]

# ╔═╡ 6719f8e0-1668-4e6e-bcfb-033e9fdf3a01
md"""
## Myocardium
"""

# ╔═╡ 7f1365a6-a73b-40a7-8b4b-bcf52eaf4161
myocardium_lacs = μ(Materials.muscle, energies)

# ╔═╡ 16ef63ae-76e2-42e0-9e9b-7e48ce4372a0
myocardium_hus = convert_hu(myocardium_lacs, energies)

# ╔═╡ 90249637-5e97-4d2b-bd80-71fbf0ac3c2d
myocardium_hus[ax]

# ╔═╡ d23b30da-3002-4e3c-b192-f6448a183e0c
md"""
## Bone
"""

# ╔═╡ b34be066-952e-4691-ac22-30832b19c527
bone_lacs = μ(Materials.corticalbone, energies)

# ╔═╡ 76d73fb9-4c86-46bc-9b45-c6ec3c19af35
bone_hus = convert_hu(bone_lacs, energies)

# ╔═╡ 54a60cb3-8b8a-44a4-937a-f4a7b879f22e
bone_hus[ax]

# ╔═╡ 3cff6264-b496-4fd6-ba40-ae9b7d407b40
md"""
## Soft tissue
"""

# ╔═╡ 07000411-c860-47e5-8331-a52d3c35e1a2
soft_tissue_lacs = μ(Materials.softtissue, energies)

# ╔═╡ 0ba53139-a282-43d6-aae4-f4e6c89d9a38
soft_tissue_hus = convert_hu(soft_tissue_lacs, energies)

# ╔═╡ d74f78c9-4928-4688-b1e7-43e34e9862f0
soft_tissue_hus[ax]

# ╔═╡ e10b6957-8a9b-477c-bb01-d193f7916ab8
md"""
## Calcium Inserts
"""

# ╔═╡ 9909e0e4-1a78-47df-a669-404036265db6
calcium_lacs = μ(Elements.Calcium, energies)

# ╔═╡ 75c1d3dd-d028-4864-bd18-6f9eadfe2fbb
const MYOCARDIUM_DENSITY = 1.050g/cm^3

# ╔═╡ 1b6ecad9-e437-4bd3-a4eb-7a12b47fef76
function calcium_mixture_lac(calcium_concentration)
	volume_fraction_calcium = calcium_concentration / MYOCARDIUM_DENSITY
	volume_fraction_myocardium = 1 - volume_fraction_calcium
	lac_mixture = volume_fraction_calcium * calcium_lacs + volume_fraction_myocardium * myocardium_lacs
	AxisArray(lac_mixture, axs)
end

# ╔═╡ ece618a5-e113-4e0a-adb0-d7728aa38c1a
begin
	insert_200_lacs = calcium_mixture_lac(0.200g/cm^3)
	insert_400_lacs = calcium_mixture_lac(0.400g/cm^3)
	insert_800_lacs = calcium_mixture_lac(0.800g/cm^3)
end

# ╔═╡ 2e4d535c-8a7b-4e1e-b647-cf1814bb1a71
begin
	insert_200_hus = convert_hu(insert_200_lacs, energies)
	insert_400_hus = convert_hu(insert_400_lacs, energies)
	insert_800_hus = convert_hu(insert_800_lacs, energies)
end

# ╔═╡ 3e6c3fd7-9831-47d6-8d0a-8044639373e8
insert_800_hus[ax]

# ╔═╡ f88ecd50-a99e-481a-8b01-888d097d29c5
md"""
# Prepare QRM Phantom
"""

# ╔═╡ 58b746d1-abce-48c7-b991-3de8db3dafc6
md"""
## Thorax
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
## Lungs
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
## Heart
"""

# ╔═╡ e81b69cc-e86f-408b-93e0-a77e2cefab1a
begin
	# Heart
	width_heart = (40mm, 40mm, 15mm)
	
	heart = cylinder(center, width_heart, angles, myocardium_hus[ax])
end;

# ╔═╡ e2819635-c8f8-48ad-a599-f6426d90dffb
md"""
## Spine
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
## Coronary Artery Calcium Inserts
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

# ╔═╡ 67edb5cf-dda7-4cac-888c-f0c1648a8612
dims = (300, 300, 20) # odd

# ╔═╡ 3a322442-93bb-447d-a672-831a075d4272
deltas = (1mm, 1mm, 1mm);

# ╔═╡ de15c717-7d9b-4b47-801a-7e8cf0352cdb
ig1 = ImageGeom(;dims = dims, deltas = deltas);

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

# ╔═╡ 8b0a6f89-4305-4e55-8f22-49d23629be95
md"""
# CT
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
## GE Lightspeed System Geometry
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
ct_geom_plot3(CtFanArc(;ps_ge_lightspeed...), ig1)

# ╔═╡ 787d1654-4c19-4d0c-9ee8-06dc461b7a95
scanner = CtFanArc(;ps_ge_lightspeed...);

# ╔═╡ 9c6343a7-2d5e-4548-ae16-8fd87af07d03
md"""
## Get Projection Data
Using `rays` from Sinograms.jl and `radon` from ImagePhantoms.jl
"""

# ╔═╡ 1eda2628-557a-4b26-ac5a-a5e1448f7518
projection_data = radon(rays(scanner), objects);

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
# ╠═1cd5f63c-36af-4c72-bccd-5c8f25dff3fe
# ╠═ae2e029b-0590-4ed4-911d-9d4a80dea2af
# ╠═fa5406de-baa2-4f66-b883-2c551ec813b8
# ╠═bb1f98c9-7163-4eaa-a1a5-8e6009a53870
# ╠═6deb0f71-a009-4ac0-b48a-ec334280eaec
# ╟─1e62720c-e9bb-4eec-a7b9-e39e6731c574
# ╠═1fc443ff-39a6-4f82-a582-6aa35ddd064f
# ╠═2d3c740d-da12-4f61-899e-ee5252e33b48
# ╠═8728f157-87c1-4a5e-a376-9f71f9e3d82a
# ╟─1f5b5c39-a125-45bd-b046-6f11060c2d78
# ╠═6e3577ae-3ea6-4ce0-9f2d-381111aa05d2
# ╠═9e7cae05-ac89-4ba7-a5c9-2a8e102ecb39
# ╠═482c2bf8-1234-4e08-ac2f-7f40dda8a1f4
# ╟─5a75e9ae-bd8e-4a15-9690-c00d95749f41
# ╠═1ed21aef-14bc-49a9-a397-c6c5d15e7c31
# ╠═8c3dcf4f-c3d5-444e-8693-e487e1ed54e5
# ╠═cc2b185d-e8b4-4283-a9d2-5af5755606d0
# ╠═9a54542a-d333-4b08-ad8e-17be90b8b0fc
# ╟─6719f8e0-1668-4e6e-bcfb-033e9fdf3a01
# ╠═7f1365a6-a73b-40a7-8b4b-bcf52eaf4161
# ╠═16ef63ae-76e2-42e0-9e9b-7e48ce4372a0
# ╠═90249637-5e97-4d2b-bd80-71fbf0ac3c2d
# ╟─d23b30da-3002-4e3c-b192-f6448a183e0c
# ╠═b34be066-952e-4691-ac22-30832b19c527
# ╠═76d73fb9-4c86-46bc-9b45-c6ec3c19af35
# ╠═54a60cb3-8b8a-44a4-937a-f4a7b879f22e
# ╟─3cff6264-b496-4fd6-ba40-ae9b7d407b40
# ╠═07000411-c860-47e5-8331-a52d3c35e1a2
# ╠═0ba53139-a282-43d6-aae4-f4e6c89d9a38
# ╠═d74f78c9-4928-4688-b1e7-43e34e9862f0
# ╟─e10b6957-8a9b-477c-bb01-d193f7916ab8
# ╠═9909e0e4-1a78-47df-a669-404036265db6
# ╠═75c1d3dd-d028-4864-bd18-6f9eadfe2fbb
# ╠═1b6ecad9-e437-4bd3-a4eb-7a12b47fef76
# ╠═ece618a5-e113-4e0a-adb0-d7728aa38c1a
# ╠═2e4d535c-8a7b-4e1e-b647-cf1814bb1a71
# ╠═3e6c3fd7-9831-47d6-8d0a-8044639373e8
# ╟─f88ecd50-a99e-481a-8b01-888d097d29c5
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
# ╠═67edb5cf-dda7-4cac-888c-f0c1648a8612
# ╠═3a322442-93bb-447d-a672-831a075d4272
# ╠═de15c717-7d9b-4b47-801a-7e8cf0352cdb
# ╠═d152f444-44b3-4f03-baed-8ea67fbd77a7
# ╠═61c53fb9-867c-48b2-91f6-b32ed721d222
# ╟─49a83445-c2f2-4d43-8dcd-343cc73fd88a
# ╟─55d7b8e6-6330-4039-a617-238a835b0daa
# ╟─8b0a6f89-4305-4e55-8f22-49d23629be95
# ╟─a8ab1a0b-7729-43b9-9044-d5238846cdf5
# ╠═14891577-24a3-4507-9e42-1c0c723e11f0
# ╟─432d5f93-70c1-4985-9f19-b066f4dc77a0
# ╠═61e34bbc-232f-499d-9343-ea944e76a760
# ╠═e3ca9f1b-1543-4b92-804b-b684170cdb70
# ╠═787d1654-4c19-4d0c-9ee8-06dc461b7a95
# ╟─9c6343a7-2d5e-4548-ae16-8fd87af07d03
# ╠═1eda2628-557a-4b26-ac5a-a5e1448f7518
# ╟─fb638c1d-e734-4b21-b32b-0bc7880355e2
# ╠═6d7ca293-2f35-418e-9211-069c9e686cb3
# ╟─dc86e054-9f05-4db6-aee3-69e379fdc26c
# ╟─5c949f69-f79d-4cf0-9dd0-698ad5dd94b2
