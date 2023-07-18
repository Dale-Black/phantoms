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

# ╔═╡ f8dd85bd-6807-4c6a-a48e-cdde9c0ef307
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

# ╔═╡ 2b60cd4a-8c2e-4f3d-855d-addf9e2fe8f2
TableOfContents()

# ╔═╡ e4773e40-abf1-4258-a6a1-0cc22b5197e0
md"""
# Prepare CT System
"""

# ╔═╡ fa588cb8-f5e3-4608-8205-77e86d260391
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

# ╔═╡ 530b2399-b6fc-4358-a191-0e840cd1d75a
md"""
## GE Lightspeed System Geometry
doi: 10.1109/TMI.2006.882141
"""

# ╔═╡ c8f96410-a7ff-479b-9526-da9b8dbd3947
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

# ╔═╡ cb4d8563-4f20-4739-8a47-1852d47ce161
ct_geom = CtFanArc(;ps_ge_lightspeed...);

# ╔═╡ f99a72c1-bb92-49c0-9ce5-fc554240241d
md"""
# Prepare Phantom
"""

# ╔═╡ 05f30df3-0ed5-47ae-96c2-5f8b46402105
md"""
## Monoenergetic Materials

The linear attenuation coefficient (often symbolized by the Greek letter μ) of a material describes how much a beam of x-ray or gamma radiation is attenuated (i.e., reduced in intensity) by that material. It is a measure of the probability per unit path length that a photon will interact with the material and lose energy.

The unit of the linear attenuation coefficient is typically inverse length (1/length), such as inverse centimeters (cm⁻¹). A higher value of the linear attenuation coefficient means the material is more effective at attenuating the radiation.

In the context of medical imaging, for example, materials like bone that have high linear attenuation coefficients appear white on an x-ray image because they absorb more of the x-ray photons, while materials like air and soft tissues that have lower coefficients appear darker because more of the x-rays pass through them.

The linear attenuation coefficient depends on the energy of the radiation, as well as the atomic number and density of the material.
"""

# ╔═╡ 97d9cae7-7577-4c32-a950-1b4893ab85e8
energies = [80keV, 100keV, 120keV, 135keV]

# ╔═╡ 5d57c14a-1f3e-4bc1-8821-7906b8eb13d1
axs = AxisArrays.Axis{:energy}(energies)

# ╔═╡ 89b42446-08aa-4220-ba97-96d09640dbb4
ax = AxisArrays.Axis{:energy}(100keV)

# ╔═╡ ceb86de9-ee21-44b0-811e-e8b025ed89a2
md"""
### Water
"""

# ╔═╡ 0ae3557e-2ae7-4199-8052-984921e9d8d5
water_lacs = μ(Materials.water, energies)

# ╔═╡ acd9421c-b3ed-4ce4-9300-4069a407a635
function convert_hu(lacs, energies)
    values = [1000 .* ((lacs[AxisArrays.Axis{:energy}(e)] ./ water_lacs[AxisArrays.Axis{:energy}(e)]) .- 1) for e in energies]
    AxisArray(values, axs)
end

# ╔═╡ 4b026893-63aa-430f-99b0-e4ee8e46d26b
water_hus = convert_hu(water_lacs, energies)

# ╔═╡ abb2a261-5275-47a1-81b2-4907b870f03b
water_hus[ax]

# ╔═╡ 660fdf4d-acc7-4724-96a2-5ba023d4c2ef
md"""
### Air
"""

# ╔═╡ e8b999bc-0620-4d14-bb16-95650c9f9903
air_lacs = μ(Materials.air, energies)

# ╔═╡ c51650ed-295e-4302-876b-4e5f017fe596
air_hus = convert_hu(air_lacs, energies)

# ╔═╡ e60b4b3c-87d5-4b5f-8037-24222904d3c2
air_hus[ax]

# ╔═╡ 1e9ab450-f06e-4f95-afe4-1be4315f56a4
md"""
### Lung
"""

# ╔═╡ 466e8930-e614-4a1a-99d5-cb8349229faf
lung_tissue_lacs = μ(Materials.lung, energies)

# ╔═╡ fc976163-e618-4fbf-9aa5-a709ff4e8df4
lung_lacs = AxisArray((0.75*air_lacs) + (0.25*lung_tissue_lacs), axs)

# ╔═╡ a5de10b3-7aed-48e2-8bb2-8c3dd63e2b27
lung_hus = convert_hu(lung_lacs, energies)

# ╔═╡ 5f09246f-a78b-4bea-8f4f-025beaaaee25
lung_hus[ax]

# ╔═╡ 0b661db0-597b-4816-be93-e33ecc8363af
md"""
### Myocardium
"""

# ╔═╡ 9420d6da-a821-491a-9988-18c96d4e7fda
myocardium_lacs = μ(Materials.muscle, energies)

# ╔═╡ 9cbcb88b-bde6-4947-9b48-8949eba876a1
myocardium_hus = convert_hu(myocardium_lacs, energies)

# ╔═╡ 9aea92e4-3652-48d5-9c63-b3d9c9c240f0
myocardium_hus[ax]

# ╔═╡ 7d6d7d79-3695-47e9-9bb5-9c3b00dc1a07
md"""
### Bone
"""

# ╔═╡ fbe77787-433b-4c7d-aa5d-3025264064c1
bone_lacs = μ(Materials.corticalbone, energies)

# ╔═╡ a33b2176-d304-4f9c-8b5a-d638225b0120
bone_hus = convert_hu(bone_lacs, energies)

# ╔═╡ 6f4d5eeb-13ab-432f-b56e-de7fbb823a73
bone_hus[ax]

# ╔═╡ d7135fbb-41d6-4d19-9185-c7c61f201e1f
md"""
### Soft tissue
"""

# ╔═╡ e445750c-55f9-44af-bf8d-9e8f4df8f014
soft_tissue_lacs = μ(Materials.softtissue, energies)

# ╔═╡ bad8bd73-aa90-427c-9d5c-8630eba3e074
soft_tissue_hus = convert_hu(soft_tissue_lacs, energies)

# ╔═╡ c7c4ae19-de2a-40e4-ba9b-3919d3929278
soft_tissue_hus[ax]

# ╔═╡ 079c46b7-eb37-4a69-8ef9-b05081301ae7
md"""
### Calcium Inserts
"""

# ╔═╡ 980737a8-7ffd-4b4f-b267-24ff8a55447c
calcium_lacs = μ(Elements.Calcium, energies)

# ╔═╡ 4d31b0c2-a4e2-4ea8-9fbd-05b466350dd2
const MYOCARDIUM_DENSITY = 1.050g/cm^3

# ╔═╡ 51923cff-43ba-4dc7-9cdc-af094e0048c4
function calcium_mixture_lac(calcium_concentration)
	volume_fraction_calcium = calcium_concentration / MYOCARDIUM_DENSITY
	volume_fraction_myocardium = 1 - volume_fraction_calcium
	lac_mixture = volume_fraction_calcium * calcium_lacs + volume_fraction_myocardium * myocardium_lacs
	AxisArray(lac_mixture, axs)
end

# ╔═╡ 442b04bd-1773-4507-9fbd-cd3f98a72cd1
begin
	insert_200_lacs = calcium_mixture_lac(0.200g/cm^3)
	insert_400_lacs = calcium_mixture_lac(0.400g/cm^3)
	insert_800_lacs = calcium_mixture_lac(0.800g/cm^3)
end

# ╔═╡ 5c3e559b-9c92-485b-881f-43c430624607
begin
	insert_200_hus = convert_hu(insert_200_lacs, energies)
	insert_400_hus = convert_hu(insert_400_lacs, energies)
	insert_800_hus = convert_hu(insert_800_lacs, energies)
end

# ╔═╡ feb1e997-6733-43d6-8ca2-b950d444ab93
insert_800_hus[ax]

# ╔═╡ ffb0457d-f66a-4154-9b7c-4fcc2f754e0a
md"""
## QRM Geometry & Materials
"""

# ╔═╡ 4d188382-18ea-490e-81f2-781b70077adf
md"""
### Thorax
"""

# ╔═╡ 2bebae16-2f20-4d8b-a808-88e002de5155
begin
	# Thorax
	center = (0mm, 0mm, 0mm)
	width = (150mm, 100mm, 15mm) # x radius, y radius, height
	angles = (0, 0, 0)
	
	thorax = cylinder(center, width, angles, soft_tissue_hus[ax])
end;

# ╔═╡ 288bf3a4-541e-4989-b60c-3c876cee5f33
md"""
### Lungs
"""

# ╔═╡ 0b573d8a-4a58-43c1-be5b-63da604a1fe4
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

# ╔═╡ d01d3eab-bba2-4738-a44a-46742525b1c5
md"""
### Heart
"""

# ╔═╡ 74a3c2ee-6ce0-43a6-8c0a-f3280cbb3e8a
begin
	# Heart
	width_heart = (40mm, 40mm, 15mm)
	
	heart = cylinder(center, width_heart, angles, myocardium_hus[ax])
end;

# ╔═╡ 3ab892f0-ba98-41cd-8ee7-22027e6a6e54
md"""
### Spine
"""

# ╔═╡ 77882031-00e2-4c6d-9bf3-1ae7e9beb0eb
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

# ╔═╡ 8bce1570-1b31-43f4-a452-8d801e1a4448
md"""
### Coronary Artery Calcium Inserts
"""

# ╔═╡ 2579e4a3-2d23-4723-b87a-4a52f3a775fe
begin
	# Small Calcium Inserts (1 mm)
	width_small_insert = (0.5mm, 0.5mm, 1mm)
	
	small_insert_low_density = cylinder((5mm, 5mm, 0mm), width_small_insert, angles, insert_200_hus[ax])
	small_insert_medium_density = cylinder((-5mm, 5mm, 0mm), width_small_insert, angles, insert_400_hus[ax])
	small_insert_high_density = cylinder((0mm, -5mm, 0mm), width_small_insert, angles, insert_800_hus[ax])
end;

# ╔═╡ 32bf6b5a-a1c9-4eda-b32c-d7e9907e8f80
begin
	# Medium Calcium Inserts (1 mm)
	width_medium_insert = (1.5mm, 1.5mm, 3mm)
	
	medium_insert_low_density = cylinder((10mm, 10mm, 0mm), width_medium_insert, angles, insert_200_hus[ax])
	medium_insert_medium_density = cylinder((-10mm, 10mm, 0mm), width_medium_insert, angles, insert_400_hus[ax])
	medium_insert_high_density = cylinder((0mm, -10mm, 0mm), width_medium_insert, angles, insert_800_hus[ax])
end;

# ╔═╡ b3b197b1-6dca-4351-907e-f65bcd728146
begin
	# Large Calcium Inserts (1 mm)
	width_large_insert = (2.5mm, 2.5mm, 5mm)
	
	large_insert_low_density = cylinder((15mm, 15mm, 0mm), width_large_insert, angles, insert_200_hus[ax])
	large_insert_medium_density = cylinder((-15mm, 15mm, 0mm), width_large_insert, angles, insert_400_hus[ax])
	large_insert_high_density = cylinder((0mm, -15mm, 0mm), width_large_insert, angles, insert_800_hus[ax])
end;

# ╔═╡ c4204486-6db9-4305-a33c-ee215dcba407
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

# ╔═╡ 3a291073-5941-4ab9-b9f8-0874b7894532
md"""
# Projection, Groundtruth, & Reconstruction
"""

# ╔═╡ deba906f-1b2e-4837-b6ac-b21ca3406a32
md"""
## Get Ground Truth Phantom
"""

# ╔═╡ 24096cac-8b07-4c4e-8a99-3fd2dcffd5fe
function generate_axes(width, deltas; offset = 5mm)
	x = ustrip.(width[1] + offset)
	Δx = ustrip.(deltas[1])

	y = ustrip.(width[2] + offset)
	Δy = ustrip.(deltas[2])

	z = ustrip.(width[3] + offset)
	Δz = ustrip.(deltas[3])

	
	x_axes = (-x:Δx:x)mm
	y_axes = (-y:Δy:y)mm
	z_axes = (-z:Δz:z)mm
	return (x_axes, y_axes, z_axes)
end

# ╔═╡ 7be7dd84-6c88-4a8b-a1d7-860091876f58
groundtruth_axes = generate_axes(width, (0.5mm, 0.5mm, 0.25mm))

# ╔═╡ aeb12303-3b1a-4269-b206-bd639000aff2
groundtruth_phantom = phantom(groundtruth_axes..., objects);

# ╔═╡ 81d5f945-cdde-4219-907e-ca1b32a28384
@bind z1 PlutoUI.Slider(Base.axes(groundtruth_phantom, 3); show_value = true, default = div(size(groundtruth_phantom, 3), 2))

# ╔═╡ 525dea64-8765-4e8c-8b9c-8ac417a54f4c
let
	f = Figure()

	ax = CairoMakie.Axis(
		f[1, 1],
		title="Ground Truth Phantom"
	)
	hm = CairoMakie.heatmap!(ustrip.(groundtruth_phantom[:, :, z1]); colormap=:grays)
	Colorbar(f[:, end+1], hm)
	hidedecorations!(ax, ticks = false, ticklabels = false)
	
	f
end

# ╔═╡ f813fe7a-55c1-45ea-87d2-1d69460030d1
md"""
## Get Projection Data
Using `rays` from Sinograms.jl and `radon` from ImagePhantoms.jl

[Radon Transform](https://en.wikipedia.org/wiki/Radon_transform)
"""

# ╔═╡ 29fada9e-1dcc-4a49-b56d-30ed637c6464
projection_data = radon(rays(ct_geom), objects);

# ╔═╡ 348dc16b-6bed-4f8d-972b-3e7fe5fa0805
md"""
## Get Reconstructed Phantom
"""

# ╔═╡ 2274141c-4a54-4095-9bd2-6268a52c2521
fov = (256mm, 256mm, 160mm) # fov (change this, not dims)

# ╔═╡ b36be38b-0a81-4263-b4c6-39c06ac2ab09
deltas = (0.5mm, 0.5mm, 0.5mm); # voxel size

# ╔═╡ f7183d3c-73c4-4b97-9e6f-abf6fe43f2a9
dims = Int.(fov ./ deltas) # number of voxels

# ╔═╡ b4165662-23c3-4cae-8b9f-2ddef9d65419
im_geom = ImageGeom(;dims = dims, deltas = deltas);

# ╔═╡ c1e35ea8-3e71-48dd-bd9e-e849da794bb6
ct_geom_plot3(CtFanArc(;ps_ge_lightspeed...), im_geom)

# ╔═╡ 7fedb79d-4c32-4bce-976d-7a54eea1e114
begin
	plan = plan_fbp(ct_geom, im_geom)
	image_data = fdk(plan, projection_data)
end;

# ╔═╡ 53512fb1-2cb3-4c7b-a650-34b1a881718d
@bind z2 PlutoUI.Slider(Base.axes(image_data, 3), show_value=true, default = div(size(image_data, 3), 2))

# ╔═╡ c4677227-f4d8-45cb-9872-992dbfaebe80
let
	f = Figure()
	
	ax = CairoMakie.Axis(
		f[1, 1],
		title="Reconstructed Image Data"
	)
	CairoMakie.heatmap!(ustrip.(image_data[:, :, z2]); colormap=:grays)
	
	f
end

# ╔═╡ fcb46905-b4e9-40a1-804b-125649664db2
groundtruth_phantom[308:312, 208:212, 80]

# ╔═╡ aed5db31-dd3d-43d1-8877-1ad929ded06d
image_data[254:258, 254:258, 160]

# ╔═╡ Cell order:
# ╠═f8dd85bd-6807-4c6a-a48e-cdde9c0ef307
# ╠═2b60cd4a-8c2e-4f3d-855d-addf9e2fe8f2
# ╟─e4773e40-abf1-4258-a6a1-0cc22b5197e0
# ╟─fa588cb8-f5e3-4608-8205-77e86d260391
# ╟─530b2399-b6fc-4358-a191-0e840cd1d75a
# ╠═c8f96410-a7ff-479b-9526-da9b8dbd3947
# ╠═c1e35ea8-3e71-48dd-bd9e-e849da794bb6
# ╠═cb4d8563-4f20-4739-8a47-1852d47ce161
# ╟─f99a72c1-bb92-49c0-9ce5-fc554240241d
# ╟─05f30df3-0ed5-47ae-96c2-5f8b46402105
# ╠═97d9cae7-7577-4c32-a950-1b4893ab85e8
# ╠═5d57c14a-1f3e-4bc1-8821-7906b8eb13d1
# ╠═89b42446-08aa-4220-ba97-96d09640dbb4
# ╠═acd9421c-b3ed-4ce4-9300-4069a407a635
# ╟─ceb86de9-ee21-44b0-811e-e8b025ed89a2
# ╠═0ae3557e-2ae7-4199-8052-984921e9d8d5
# ╠═4b026893-63aa-430f-99b0-e4ee8e46d26b
# ╠═abb2a261-5275-47a1-81b2-4907b870f03b
# ╟─660fdf4d-acc7-4724-96a2-5ba023d4c2ef
# ╠═e8b999bc-0620-4d14-bb16-95650c9f9903
# ╠═c51650ed-295e-4302-876b-4e5f017fe596
# ╠═e60b4b3c-87d5-4b5f-8037-24222904d3c2
# ╟─1e9ab450-f06e-4f95-afe4-1be4315f56a4
# ╠═466e8930-e614-4a1a-99d5-cb8349229faf
# ╠═fc976163-e618-4fbf-9aa5-a709ff4e8df4
# ╠═a5de10b3-7aed-48e2-8bb2-8c3dd63e2b27
# ╠═5f09246f-a78b-4bea-8f4f-025beaaaee25
# ╟─0b661db0-597b-4816-be93-e33ecc8363af
# ╠═9420d6da-a821-491a-9988-18c96d4e7fda
# ╠═9cbcb88b-bde6-4947-9b48-8949eba876a1
# ╠═9aea92e4-3652-48d5-9c63-b3d9c9c240f0
# ╟─7d6d7d79-3695-47e9-9bb5-9c3b00dc1a07
# ╠═fbe77787-433b-4c7d-aa5d-3025264064c1
# ╠═a33b2176-d304-4f9c-8b5a-d638225b0120
# ╠═6f4d5eeb-13ab-432f-b56e-de7fbb823a73
# ╟─d7135fbb-41d6-4d19-9185-c7c61f201e1f
# ╠═e445750c-55f9-44af-bf8d-9e8f4df8f014
# ╠═bad8bd73-aa90-427c-9d5c-8630eba3e074
# ╠═c7c4ae19-de2a-40e4-ba9b-3919d3929278
# ╟─079c46b7-eb37-4a69-8ef9-b05081301ae7
# ╠═980737a8-7ffd-4b4f-b267-24ff8a55447c
# ╠═4d31b0c2-a4e2-4ea8-9fbd-05b466350dd2
# ╠═51923cff-43ba-4dc7-9cdc-af094e0048c4
# ╠═442b04bd-1773-4507-9fbd-cd3f98a72cd1
# ╠═5c3e559b-9c92-485b-881f-43c430624607
# ╠═feb1e997-6733-43d6-8ca2-b950d444ab93
# ╟─ffb0457d-f66a-4154-9b7c-4fcc2f754e0a
# ╟─4d188382-18ea-490e-81f2-781b70077adf
# ╠═2bebae16-2f20-4d8b-a808-88e002de5155
# ╟─288bf3a4-541e-4989-b60c-3c876cee5f33
# ╠═0b573d8a-4a58-43c1-be5b-63da604a1fe4
# ╟─d01d3eab-bba2-4738-a44a-46742525b1c5
# ╠═74a3c2ee-6ce0-43a6-8c0a-f3280cbb3e8a
# ╟─3ab892f0-ba98-41cd-8ee7-22027e6a6e54
# ╠═77882031-00e2-4c6d-9bf3-1ae7e9beb0eb
# ╟─8bce1570-1b31-43f4-a452-8d801e1a4448
# ╠═2579e4a3-2d23-4723-b87a-4a52f3a775fe
# ╠═32bf6b5a-a1c9-4eda-b32c-d7e9907e8f80
# ╠═b3b197b1-6dca-4351-907e-f65bcd728146
# ╠═c4204486-6db9-4305-a33c-ee215dcba407
# ╟─3a291073-5941-4ab9-b9f8-0874b7894532
# ╟─deba906f-1b2e-4837-b6ac-b21ca3406a32
# ╠═24096cac-8b07-4c4e-8a99-3fd2dcffd5fe
# ╠═7be7dd84-6c88-4a8b-a1d7-860091876f58
# ╠═aeb12303-3b1a-4269-b206-bd639000aff2
# ╟─81d5f945-cdde-4219-907e-ca1b32a28384
# ╟─525dea64-8765-4e8c-8b9c-8ac417a54f4c
# ╟─f813fe7a-55c1-45ea-87d2-1d69460030d1
# ╠═29fada9e-1dcc-4a49-b56d-30ed637c6464
# ╟─348dc16b-6bed-4f8d-972b-3e7fe5fa0805
# ╠═2274141c-4a54-4095-9bd2-6268a52c2521
# ╠═b36be38b-0a81-4263-b4c6-39c06ac2ab09
# ╠═f7183d3c-73c4-4b97-9e6f-abf6fe43f2a9
# ╠═b4165662-23c3-4cae-8b9f-2ddef9d65419
# ╠═7fedb79d-4c32-4bce-976d-7a54eea1e114
# ╟─53512fb1-2cb3-4c7b-a650-34b1a881718d
# ╟─c4677227-f4d8-45cb-9872-992dbfaebe80
# ╠═fcb46905-b4e9-40a1-804b-125649664db2
# ╠═aed5db31-dd3d-43d1-8877-1ad929ded06d
