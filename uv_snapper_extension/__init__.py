import bpy
import bmesh
import numpy
import copy
import math
from timeit import default_timer as timer

bpy.types.Scene.atlas = bpy.props.PointerProperty(type=bpy.types.Object)


#Custom data class to hold info about the atlas box that the UV selection is in
class Box:
	def __init__(self, bounds):
		self.left_bound = bounds[0]
		self.top_bound = bounds[1]
		self.right_bound = bounds[2]
		self.bottom_bound = bounds[3]
		self.tl_corner = numpy.array((self.left_bound, self.top_bound))
		self.tr_corner = numpy.array((self.right_bound, self.top_bound))
		self.br_corner = numpy.array((self.right_bound, self.bottom_bound))
		self.bl_corner = numpy.array((self.left_bound, self.bottom_bound))
		self.centre = numpy.array((((self.right_bound + self.left_bound)/2), ((self.top_bound + self.bottom_bound)/2)))
		self.length = abs(self.top_bound - self.bottom_bound)

PROPS = [ 
	('use_full_soft_select', bpy.props.BoolProperty(name='Use Full Soft Select', default = True )),
	('padding', bpy.props.FloatProperty(name = 'Shell Padding Divisor', default = 50)),
	('softselect_falloff', bpy.props.FloatProperty(name = 'Soft Select Falloff Range', default = 1, soft_max = 1, min = 0)),
	('softselect_multiplier', bpy.props.FloatProperty(name = 'Soft Select Multiplier', default = 1)),
	]

#UI Class
class UV_Snapper_UV_PT_Panel(bpy.types.Panel):
	bl_idname = "UV_Snapper_UV_PT_Panel"
	bl_label = "UV Snapper"
	bl_space_type = "IMAGE_EDITOR"
	bl_region_type = "UI"
	bl_category = "UV Snapper"


	def draw(self, context):

		layout = self.layout
		box = layout.box()
		col = box.column(align = True)
		row = col.row (align =True)

		obj = bpy.context.active_object

		col.operator("uv.snaptoatlas", text = "Snap UV to Atlas")
		if obj.mode == "EDIT":

			col.enabled = True
		else:
			col.enabled = False

		col = box.column(align = True)
		col.prop_search(context.scene, "atlas", context.scene, "objects", text="Atlas")

			 

class UV_Snapper_Settings_UV_PT_Panel(bpy.types.Panel):
	bl_idname = "UV_Snapper_Settings_UV_PT_Panel"
	bl_label = "Settings"
	bl_space_type = "IMAGE_EDITOR"
	bl_region_type = "UI"
	bl_parent_id = "UV_Snapper_UV_PT_Panel"
	bl_category = "UV Snapper"
	bl_options = {'DEFAULT_CLOSED'}

	def draw(self, context):

		layout = self.layout
		box = layout.box()
		col = box.column(align = True)
		row = col.row (align =True)


		col = box.column(align = True)
		col.prop(context.scene, "use_full_soft_select")
		col = box.column(align = True)
		col.prop(context.scene, "padding")
		col = box.column(align = True)
		col.prop(context.scene, "softselect_falloff")
		col = box.column(align = True)
		col.prop(context.scene, "softselect_multiplier")
			 	   



#---GENERAL USE FUNCTIONS----

def FindBoundsofVector2Array(uv_array): #takes an array of V2s and finds the highest and lowest values in each axis

	default_uv = uv_array[0]
	bounds_array = [default_uv[0],default_uv[1],default_uv[0],default_uv[1]] #output is low u, high v, high u, low v (clockwise from left)

	bounds_array[0] = min(uv[0] for uv in uv_array)	
	bounds_array[1] = max(uv[1] for uv in uv_array)
	bounds_array[2] = max(uv[0] for uv in uv_array)
	bounds_array[3] = min(uv[1] for uv in uv_array)

	return bounds_array

def FindClosestVertfromArray( single_vector, vert_array, rounding, numpy): #takes a single v2 and an array of verts, finds the closest vert to the point, returns vert

	vert_dist_array = []

	for vert in vert_array:	#loop through all verts, find the distance to point
		loc = vert.co.xy
		dist = numpy.linalg.norm(single_vector - loc)

		#this controls if we round up, down or in both directions. in any direction we only take positive values, if we round down we flip signs
		if rounding == 0:		
			dist = abs(dist)
			vert_dist_array.append((vert, dist))
		if rounding == 1:		
			if dist >= 0:
				vert_dist_array.append((vert, dist))
		if rounding == -1:
			dist = dist * -1
			if dist >= 0:
				vert_dist_array.append((vert, dist))

	output_vert = min(vert_dist_array, key=lambda i: i[1])[0]	#gets the vert with the lowest distance

	return output_vert

def FindClosestValuefromArray(single_value, values_array, rounding): #same as FindClosestVert but for float values 

	value_dist_array = []

	for value in values_array:
		dist = value - single_value 
		if rounding == 0:		
			dist = abs(dist)
			value_dist_array.append((value, dist))
		if rounding == 1:		
			if dist >= 0:
				value_dist_array.append((value, dist))
		if rounding == -1:
			dist = dist * -1
			if dist >= 0:
				value_dist_array.append((value, dist))

	closest_value = min(value_dist_array, key=lambda i: i[1])[0]

	return closest_value

def FindCentreofBounds(bounds_array):	#finds the centre point of a bounds array

	averageU = (bounds_array[0] + bounds_array[2])/2	
	averageV = (bounds_array[1] + bounds_array[3])/2

	centre = (averageU, averageV)

	return centre 

def Lerp(a, b, t):
	return (a +(b-a)*t) 

def InverseLerp(a,b ,t):
	return (t-a)/(b-a) 

def GetSampledBoundsPosition(bound, vert_pos, axis):  # returns a value from the opposite axis provided based on provided bounds 

	#if the left x bound is 0,0 at the bottom and 0.5, 1 at the top, and we feed in a y value of 0.5 we get back 0.25 (the min x value it should have at that y value)

	if axis == 0:	#gives the point on the left and right bounds (x coords) that corresponds to the Y coordinate of the vert
		point_on_axis = InverseLerp(bound[0][1], bound[1][1], vert_pos[1])
		sampled_bounds_pos = Lerp(bound[0][0], bound[1][0], point_on_axis)
	else:			#gives the point on the top and bottom bounds that corresponds to the X coordinate of the vert
		point_on_axis = InverseLerp(bound[0][0], bound[1][0], vert_pos[0])
		sampled_bounds_pos = Lerp(bound[0][1], bound[1][1], point_on_axis)

	return sampled_bounds_pos 

def OutOfBounds(coords, box):	#checks if a point is outwith defined bounds
	bound_value = 0
	bound_direction = 0
	distance_oob = 0
	bound_axis = None

	if coords[0] > box.right_bound:
		distance_oob = coords[0] - box.right_bound
		bound_value = box.right_bound
		bound_direction = 1
		bound_axis ="X"
	if coords[0] < box.left_bound:
		if abs(coords[0] - box.left_bound) > abs(distance_oob):
			distance_oob = coords[0] - box.left_bound
			bound_value = box.left_bound
			bound_direction = -1
			bound_axis ="X"
	if coords[1] > box.top_bound:
		if abs(coords[1] - box.top_bound) > abs(distance_oob):
			distance_oob = coords[1] - box.top_bound
			bound_value = box.top_bound
			bound_direction = 1
			bound_axis ="Y"
	if coords[1] < box.bottom_bound:
		if abs(coords[1] - box.bottom_bound) > abs(distance_oob):
			distance_oob = coords[1] - box.bottom_bound
			bound_value = box.bottom_bound
			bound_direction = -1
			bound_axis ="Y"

	return distance_oob, bound_value, bound_direction, bound_axis #determines if a point is out of bounds of the current box, and if so by how much etc

def Saturate(value): #0-1 clamp
	if value < 0:
		return 0
	if value > 1:
		return 1 
	else:
		return value

def ScaleShell (bm, UVBounds, vert_array, CurrentBox, context): #scales the shell so that it fills the atlas box, plus a defined margin


	xlocs_array = []
	ylocs_array = []
	
	uvcentre = FindCentreofBounds(UVBounds)
	for vert in vert_array:
		xlocs_array.append(vert.co.x)
		ylocs_array.append(vert.co.y)

	lowest_x = (xlocs_array.index(min(xlocs_array)), min(xlocs_array))
	highest_x = (xlocs_array.index(max(xlocs_array)), max(xlocs_array))
	lowest_y = (ylocs_array.index(min(ylocs_array)), min(ylocs_array))
	highest_y = (ylocs_array.index(max(ylocs_array)), max(ylocs_array))

	dist_low_x = abs(CurrentBox.centre[0] - lowest_x[1])
	dist_high_x = abs(CurrentBox.centre[0] - highest_x[1])

	dist_low_y = abs(CurrentBox.centre[1] - lowest_y[1])
	dist_high_y = abs(CurrentBox.centre[1] - highest_y[1])

	padding = CurrentBox.length/context.scene.padding #user defined variable, defaults to 50

	if dist_low_x > dist_high_x:

		furthest_x_dist = dist_low_x
		furthest_x_vert = lowest_x
		xbound = CurrentBox.left_bound + padding
	else:
		furthest_x_dist = dist_high_x
		furthest_x_vert = highest_x
		xbound = CurrentBox.right_bound - padding

	if dist_low_y > dist_high_y:
		furthest_y_dist = dist_low_y
		furthest_y_vert = lowest_y
		ybound = CurrentBox.bottom_bound + padding
	else:
		furthest_y_dist = dist_high_y
		furthest_y_vert = lowest_y
		ybound = CurrentBox.top_bound - padding

	#scale to get the bounds to the desired location
	xscale = abs(xbound - CurrentBox.centre[0])/ furthest_x_dist 
	yscale = abs(ybound - CurrentBox.centre[1])/ furthest_y_dist

	for vert in vert_array:
		xdist_from_centre = vert.co.x - CurrentBox.centre[0]
		vert.co.x = CurrentBox.centre[0] + (xdist_from_centre * xscale)
		ydist_from_centre = vert.co.y - CurrentBox.centre[1]
		vert.co.y = CurrentBox.centre[1] + (ydist_from_centre * yscale)

	return bm

def CentreShell(bm, UVBounds, vert_array, CurrentBox, context): #moves the current shell to the centre of the selected box


	uvcentre = FindCentreofBounds(UVBounds)
	dist_between_centres = CurrentBox.centre - uvcentre
	for vert in vert_array:
		vert.co.x += dist_between_centres[0]
		vert.co.y += dist_between_centres[1]

	return bm

#---STAGES OF EXECUTION----

def ReadAtlas(context): #reads data from the atlas object

	#define atlas object to draw loops from
	atlas_mesh = context.scene.atlas.data 

	#create a bmesh from the atlas_mesh to grab uvs
	bm = bmesh.new()
	bm.from_mesh(atlas_mesh)						
	uvlayer = bm.loops.layers.uv.verify()
								
	atlas_u_coords = []
	atlas_v_coords = []

	#iterate through every uv, and append it's coords to the list if not already there
	for face in bm.faces:								
		for loop in face.loops:
			if loop[uvlayer].uv[0] not in atlas_u_coords:
				atlas_u_coords.append(loop[uvlayer].uv[0])
			if loop[uvlayer].uv[1] not in atlas_v_coords:
				atlas_v_coords.append(loop[uvlayer].uv[1])

	return atlas_u_coords, atlas_v_coords

def SetupBmesh(context): #splits off the originally selected uvs on the actual mesh, sets up a bmesh copy of the original for all further operations


	#we need uv sync off immediately, so we save it and toggle here
	og_uv_sync = bpy.context.tool_settings.use_uv_select_sync
	bpy.context.tool_settings.use_uv_select_sync = False

	if og_uv_sync == True: #when turning sync off we usually lose the whole selection, so select all grabs everything that was selected
		bpy.ops.uv.select_all(action='SELECT')

	og_obj = bpy.context.active_object

	#ensure we're in edit mode
	bpy.ops.object.mode_set(mode = 'EDIT')

	#create a temporary bmesh of the original object
	temp_bm = bmesh.from_edit_mesh(bpy.context.active_object.data)
	temp_uvlayer = temp_bm.loops.layers.uv.verify()

	#valdiate uv selection
	temp_bm, selected_uv_coords, selected_uv_loops = FindSelectedUVs(temp_bm, temp_uvlayer, context)


	#for all faces with all loops selected, move those loops by a tiny offset to separate them from their stacked siblings
	for loop in selected_uv_loops:
		loop[temp_uvlayer].uv.x += 0.0001
		loop.tag = True


	#update bmesh applies the transforms so we can discard the bmesh 
	bmesh.update_edit_mesh(bpy.context.active_object.data)
	temp_bm.free()

	#duplicate the original mesh (not bmesh!) and do all operations to the duplicate to protect the original
	bpy.ops.object.mode_set(mode = 'OBJECT')
	bpy.ops.object.duplicate()                #using the duplicate operator because it automatically selects the duplicate
	bpy.ops.object.editmode_toggle()
	bpy.context.view_layer.objects.active.name = "uvsnapper_temp"

	#set up the duplicate as a bmesh
	bm = bmesh.from_edit_mesh(bpy.context.active_object.data)						
	uvlayer = bm.loops.layers.uv.verify()
	

	return bm, uvlayer, og_obj, og_uv_sync

def SaveInitialState(bm, uvlayer): #saves some states that will be changed later so we can reinstate them at the end

	#keeps the mesh select mode (face/edge/vert) so it can be reinstated later 
	og_mesh_select_mode = bpy.context.tool_settings.mesh_select_mode[:]

	#saves preexisting seams to be reinstated, removes them (we use seam selection for something else later)
	og_uv_seams = []
	for edge in bm.edges:
		if edge.seam == True:
			og_uv_seams.append(edge.index)
			edge.seam = False


	#grabs the indices of the originally selected mesh elements
	og_selection = [[],[],[]]

	if og_mesh_select_mode[0] == True:
		for vert in bm.verts:
			if vert.select == True:
				og_selection[0].append(vert.index)
	if og_mesh_select_mode[1] == True:
		for edge in bm.edges:
			if edge.select == True:
				og_selection[1].append(edge.index)
	if og_mesh_select_mode[2] == True:
		for face in bm.faces:
			if face.select == True:
				og_selection[2].append(face.index)

	og_uv_select_mode = bpy.context.tool_settings.uv_select_mode

	return og_mesh_select_mode, og_uv_seams, og_selection, og_uv_select_mode

def FindSelectedUVs(bm, uvlayer, context): #validate selected uvs (removes UVs that are not visible but still selected)  

	selected_uv_coords = []
	original_selection = []

	coords_to_pinned = {}

	#remove and store pins for all uvs
	for face in bm.faces:
		for loop in face.loops:
			coords_to_pinned[tuple(loop[uvlayer].uv)] = copy.deepcopy(loop[uvlayer].pin_uv)
			loop[uvlayer].pin_uv = False


	#go through all faces and check if they are selected
	for face in bm.faces:
		selected_face = True								
		for loop in face.loops:
			if loop[uvlayer].select == False:
				selected_face = False
				break

		#if all the loops of a face are selected then the face itself is - pin all its loops
		if selected_face == True:
			#pin all selected uvs, add them to an array
			for loop in face.loops:
				original_selection.append(loop)
				loop[uvlayer].pin_uv = True

	#select all visible, remove pin (this is used to figure out which uvs are visible since the select all operator only works on visible)
	bpy.ops.uv.select_all(action='SELECT')
	bpy.ops.uv.pin(clear=True)
	bpy.ops.uv.select_all(action='DESELECT')

	selected_uv_loops = []
	#all uvs that were in the original selection are now unpinned. Add them to the output array
	for loop in original_selection:
		if loop[uvlayer].pin_uv == False: 
			selected_uv_coords.append(copy.deepcopy(loop[uvlayer].uv))	#deep copy is needed here to prevent interefence further down the line
			selected_uv_loops.append(loop)
			loop[uvlayer].select = True
			loop.tag = True       			#all valid loops are tagged for the next stage

	for face in bm.faces:
		for loop in face.loops:
			loop[uvlayer].pin_uv = coords_to_pinned[tuple(loop[uvlayer].uv)]	#restore initial pinned state

	return bm, selected_uv_coords, selected_uv_loops

def CreateWorkingDuplicate(bm, SelectedUVs, uvlayer, og_uv_seams): #creates the duplicate mesh that we use to do uv operations in using mesh elements

	bpy.context.tool_settings.uv_select_mode = "VERTEX"
	og_uv_locs = {}		#used later to match up the translated verts with their original coords so we can copy them back at the end

	unselected_faces = []
	selected_faces = []

	for face in bm.faces:			#should be able to get this data from the FindSelectedUVs stage?
		face_good = True
		for loop in face.loops:
			if loop.tag == True:
				continue
			else:
				face_good = False
		if face_good == True:
			selected_faces.append(face)
			face.select = True
		else:
			unselected_faces.append(face)

	bmesh.ops.delete( bm, geom = unselected_faces, context = 'FACES' )

	# ----- This section is required to get the geometry to match the UV layout. It serves to split vertices that have more than one UV position

	sep_edges = []
	bpy.ops.uv.seams_from_islands() #we use this so that we don't catch random edges from looping round the side of shells
	for vert in bm.verts: #this is now a reasonable set since we delete what we don't need 
		loc_to_edges = {}
		loc_array = [] 
		   
		for loop in vert.link_loops:
			if loop[uvlayer].uv not in loc_array:		#stores unique uv coords  #having issues with it wrapping round shells
				loc_array.append(loop[uvlayer].uv)
				loc_to_edges[tuple(loop[uvlayer].uv)] = []		#adds a blank array of this uv coord in the dict
		
		if len(loc_array) > 1:					#if this vert has more than one uv coords:
			for loop in vert.link_loops:
				if loop[uvlayer].uv in loc_array:
					face = loop.face
					for edge in face.edges:	
						loc_to_edges[tuple(loop[uvlayer].uv)].append(edge)	#all edges on the face the loop is attached to (all faces surrounding vertex?)
			
			edges_array = []
			
			for loc in loc_array:							#for each unique uv location for the vertex:
				for edge in loc_to_edges[tuple(loc)]:
					if edge not in edges_array:
						edges_array.append(edge)	
					else:
						if edge not in sep_edges and edge.seam == True:			#if an edge is found twice it should be split 
							sep_edges.append(edge)

	bmesh.ops.split_edges(bm, edges = sep_edges)

	bpy.ops.mesh.select_all(action='DESELECT')

	vert_array = [] #this is just bm.verts really. Think this used to have a separate purpose but is a bit redundant now

	if len(bm.verts) <= 4:
		print("four or fewer verts found. Shell marked as corner only")	#This allows the shell to skip some steps further on
		is_shell_only_corners = True
	else:
		is_shell_only_corners = False

	#----------------- This is the stage where we move the verts in 3D to create a geometric representation of the UVs.
	for vert in bm.verts:								
		og_coords = vert.link_loops[0][uvlayer].uv 
		og_uv_locs[vert] = og_coords 			#store the starting uv coords ina dict

		#sets the position of each vertex in 3D space to it's UV coordinates
		vert.co.x = og_coords[0] 	
		vert.co.y = og_coords[1]
		vert.co.z = 0	
		vert_array.append(vert)
	
	for edge in bm.edges:
		if len(edge.link_faces) < 2: #this isolates the edge of the mesh
			edge.seam = True

	#different selection sets we need for now 
	edge_verts_array = [] 
	edge_uv_array = []

	#populate a dict with the original location of each edge vert, using the index as the key - need this for doing edge snapping after soft select
	og_edge_verts_locs = {}


	#-----------WE ARE NOW WORKING IN OBJECT SPACE, NOT UV!------------------

	#select all seam edges (boundary verts of the shell), fill in selection sets  
	for edge in bm.edges:
		if edge.seam == True:
			edge.select = True
			for vert in edge.verts:
				if vert not in edge_verts_array:
					edge_verts_array.append(vert)
					vert.select = True
					edge_uv_array.append(vert.co)
					og_edge_verts_locs[vert.index] = copy.deepcopy(vert.co)	#this is probably redundant

	#we stick in vertex selection mode for the rest of the process 
	bpy.ops.mesh.select_mode(type='VERT')	

	return bm, og_uv_locs,edge_verts_array, edge_uv_array, vert_array, og_edge_verts_locs, is_shell_only_corners

def FindCurrentBox(bm, vert_array, edge_uv_array, atlas_u, atlas_v): #finds the box that the selected uvs lie within


	#Find the bounding box of the verts
	UVBounds = FindBoundsofVector2Array(edge_uv_array)	
	UVBoundsCentre = FindCentreofBounds(UVBounds)

	#this prevents a break if the shell has already been snapped exactly to the bounds - we just scale it down slightly first
	for vert in vert_array:
		xdist_from_centre = vert.co.x - UVBoundsCentre[0]
		vert.co.x = UVBoundsCentre[0] + (xdist_from_centre * 0.99)
		ydist_from_centre = vert.co.y - UVBoundsCentre[1]
		vert.co.y = UVBoundsCentre[1] + (ydist_from_centre * 0.99)

	UVBounds = FindBoundsofVector2Array(edge_uv_array)	
	UVBoundsCentre = FindCentreofBounds(UVBounds)

	#Find which atlas box the bounds lie within - Left, Top, Right, Bottom

	boxbounds = []

	boxbounds.append(FindClosestValuefromArray(UVBounds[0], atlas_u, -1))
	boxbounds.append(FindClosestValuefromArray(UVBounds[1], atlas_v, 1))
	boxbounds.append(FindClosestValuefromArray(UVBounds[2], atlas_u, 1))
	boxbounds.append(FindClosestValuefromArray(UVBounds[3], atlas_v, -1))

	CurrentBox = Box(boxbounds) #custom data type, calculates a few useful pieces of data

	return CurrentBox, UVBoundsCentre, UVBounds 

def FindAndSnapCorners(bm, uvlayer, vert_array, CurrentBox, outline_verts, UVBounds, context):

	bm = CentreShell(bm, UVBounds, vert_array, CurrentBox, context)  #move shell into centre of bxo

	box_corner_array = [CurrentBox.tl_corner, CurrentBox.tr_corner, CurrentBox.bl_corner, CurrentBox.br_corner]	#put the box data into an array for ease of use
	pre_shellscale_corner_locs = {} #corner posiitons after the bounds correction but before snapping
	pre_correction_corner_locs = {} #corner positions before anything has been done to them
	corners = []
	proximity_array = []

	usable_verts = copy.copy(outline_verts)

	num_corners = 4

	#find corners initially to allow scale shell to work properly
	for i in range(4):
		if len(usable_verts) < 1:	#this covers 3 vert shells
			num_corners = 3
			corners.append(corners[0])	#4th corner is the same as 1st
		else:
			corners.append(FindClosestVertfromArray(box_corner_array[i], usable_verts, 0, numpy))	#find the closest vertex to each corner
			pre_correction_corner_locs[corners[i]] = corners[i].co.xy 
			usable_verts.remove(corners[i])		

			
	#scale shell up so that it fills the box better

	left_ref_bound = [corners[0], corners[2]]
	top_ref_bound = [corners[0], corners[1]]
	right_ref_bound = [corners[1], corners[3]]
	bottom_ref_bound = [corners[2], corners[3]]

	corner_bounds_reference = [left_ref_bound, top_ref_bound, right_ref_bound, bottom_ref_bound]

	if num_corners == 4:

		left_bound = [corners[0].co.xy, corners[2].co.xy]
		top_bound = [corners[0].co.xy, corners[1].co.xy]
		right_bound = [corners[1].co.xy, corners[3].co.xy]
		bottom_bound = [corners[2].co.xy, corners[3].co.xy]

		corner_bounds = [left_bound, top_bound, right_bound, bottom_bound]

		left_oob = []
		top_oob = []
		right_oob =[]
		bottom_oob = []
		oob_array = [left_oob, top_oob, right_oob, bottom_oob]  #all the verts that are out of bounds in each axis

		#move corners so that all verts lie within a bounding box (all the corners joined up). This helps with keeping soft select within the atlas box
		#This ended up being quite long due to having to hardcode axis and direction variables. Tried doing it using loops but it was more trouble than it was worth
		i = -1
		for bound in corner_bounds:
			i += 1
			for vert in usable_verts:
				if i == 0 or i == 2: 																		#left/right bounds - check for x
					sampled_bounds_pos = GetSampledBoundsPosition(bound,vert.co, 0)
					if i == 2: 																						#out of right bounds
						if vert.co.x > sampled_bounds_pos:
							right_oob.append(vert)
					else:
						if vert.co.x < sampled_bounds_pos:																#out of left bounds
							left_oob.append(vert)
				else: 																							#top/bottom bounds - check for y
					sampled_bounds_pos = GetSampledBoundsPosition(bound, vert.co, 1)
					if i == 1:
						if vert.co.y > sampled_bounds_pos:															#out of bounds top
							top_oob.append(vert)
					else:		
						if vert.co.y < sampled_bounds_pos:															#out of bounds bototm
							bottom_oob.append(vert)

		it = -1

		left_bound_offset = 0
		top_bound_offset = 0
		right_bound_offset = 0
		bottom_bound_offset = 0

		bound_offset_array = [left_bound_offset,top_bound_offset, right_bound_offset, bottom_bound_offset] #how far the most out of bound vert is in each axis

		for array in oob_array:
			it += 1
			bound = corner_bounds[it]
			if len(array)> 0:
				distance_values = []
				if it == 0 :				#left bound
					axis = 0
					direction = -1
				if it == 1:					#top bound
					axis = 1
					direction = 1
				if it == 2: 				#right bound
					axis = 0
					direction = 1
				if it == 3:					#bottom bound
					axis = 1
					direction = -1

				for vert in array:
					distance_values.append(vert.co.xy)
				if direction == 1:
					c_vert = max(distance_values, key=lambda item: item[axis])	#get highest out bounds
				else:
					c_vert = min(distance_values, key=lambda item: item[axis])	#get lowest out of bounds (negative, for left and bottom bounds)
				sampled_bounds_pos = GetSampledBoundsPosition(bound , c_vert, axis)		#get the sample position for the most out of bound vert
				bound_offset_array[it] = (c_vert[axis] - sampled_bounds_pos)	#the offset is the amount that vert is out of bounds by

		#move the corners so verts are no longer out of bounds
		i = -1
		for offset in bound_offset_array:
			i += 1
			for corner in corner_bounds_reference[i]:
				if i == 0 or i == 2:
					corner.co.x += offset
				else:
					corner.co.y += offset

	#save out the corner positions before scaling or snapping
	for i in range(num_corners):	
		pre_shellscale_corner_locs[corners[i].index] = copy.copy(corners[i].co.xy)

	corner_locs = []
	for corner in corners:
		corner_locs.append(corner.co)

	#recalculate UVBounds for any future operations
	UVBounds = 	FindBoundsofVector2Array(corner_locs)
	
	pre_snap_corners_locs= {}

	#vertex locations before shell scaling
	vert_og_locs_dict = {}
	for vert in vert_array:
		vert_og_locs_dict[vert] = vert.co.xy

	#Scale the shell so that it fits the box as best as possible - reduces the work soft select has to do
	bm = ScaleShell(bm, UVBounds, vert_array, CurrentBox, context)

	#save out corner positions after scaling but before snapping (for soft select)
	for i in range(num_corners):	
		pre_snap_corners_locs[corners[i].index] = corners[i].co.xy


	for i in range(num_corners):												
		proximity_array.append(abs(numpy.linalg.norm(numpy.array(corners[i].co.xy) - box_corner_array[i])))	#distance between corner and box corner
		corners[i].co.xy = box_corner_array[i]		#snap corner positions to box corners

	closest_ind = proximity_array.index(min(proximity_array))	 #find the closest corner vertex to a box corner so we can use that as the starting point for edge marching


	return vert_array, corners, pre_shellscale_corner_locs, closest_ind, UVBounds, pre_correction_corner_locs, vert_og_locs_dict, pre_snap_corners_locs #finds the point closest to each of the boxes corners, and snaps them to said corners

def MarchEdges(bm, corners, edge_verts_array,  atlas_u, atlas_v, UVBoundsCentre, starting_vert, og_edge_verts_locs, CurrentBox): #passes through all edge vertices and snaps them to the box edges

	#This section is quite old and I can't remember some specifics

	used_verts = []
	group0 = []
	group1 = []
	group2 = []
	group3 = []
	edge_groups = [group0, group1, group2, group3]
	corner_index = -1

	#get all edge verts that aren't corners
	snappable_verts = []
	for vert in edge_verts_array:
		if vert not in corners:
			snappable_verts.append(vert)

	#loop through all edge verts and corners
	for i in range(len(edge_verts_array) + 4): #starting_vert is in the current vert in the marcher here

		if starting_vert in corners:			#if we hit a corner then increment the corner index
			corner_index += 1  
		if corner_index > 3:					#if corner index >3 we've looped fully so stop
			continue

		if starting_vert not in corners and starting_vert not in used_verts:		#if the vert isn't a corner then we can add it to the current snapping array
			edge_groups[corner_index].append(starting_vert)

		for edge in starting_vert.link_edges:	#this is the marcher part
			if edge.seam == True and edge.other_vert(starting_vert) not in used_verts:
				next_vert = edge.other_vert(starting_vert) 

		used_verts.append(starting_vert)		#save used verts so we don't backtrack
		starting_vert = next_vert				#increment

	for vert in bm.verts:	#select only edge non corner verts, debug purposes?
		vert.select = False
		if vert in used_verts:
			vert.select = True	

	snap_direction = None		#initialise snap direction variable

	group_aspect_array = []	#used to figure out which direction is the most likely to be correct, we proceed in lockstep from there
	group_lengths_array = []
	group_centres_array = []

	for group in edge_groups:	#four groups

		#we're calculating the bounds of the edge group to find it's aspect ratio, so that we can use the highest aspect ratio as the best guess to start us off 
		if len(group) > 1:	#if there is more than one vert in the group we can calculate the bounds and centre normally	

			original_loc_array = []
			for vert in group:
				original_loc_array.append(og_edge_verts_locs[vert.index])	#location before soft select, more accurate to original selection


			bounds = FindBoundsofVector2Array(original_loc_array) #calc bounds using original locations of the verts

			centre = FindCentreofBounds(bounds)
			u_length = abs(bounds[2] - bounds[0])
			v_length = abs(bounds[1] - bounds[3])


		elif len(group) == 1:	#if the length is one (one vert between two corners) calculate bounds differently	

			vert = group[0]
			original_vert_loc = og_edge_verts_locs[vert.index]	#find the original location of the single vert 
			edge_group_og_locs = [original_vert_loc]				#make an array for the original location of the single vert and it's two adjacent corners

			for e in vert.link_edges:								#find the two adjacent corners and add original loc to array
				v = e.other_vert(vert)
				if v in corners:
					edge_group_og_locs.append(og_edge_verts_locs[v.index])

			bounds = FindBoundsofVector2Array(edge_group_og_locs)		#calculate the bounds from the three verts original locs. This ensures the snapping axis is accurate
			u_length = abs(bounds[2] - bounds[0])
			v_length = abs(bounds[3] - bounds[1])
			centre = vert.co  										#use the current position of the vert to determine which direction we're snapping (rounding up or down)

		else:														#this covers off the case of a zero length group (straight line)
			u_length = 0.001
			v_length = 0.001
			centre = (0.0001, 0.0001)

		if u_length > v_length:
			if v_length != 0:
				aspect_ratio = (u_length/v_length)
			else:
				aspect_ratio = 1000
		else: 
			if u_length != 0:
				aspect_ratio = (v_length/u_length)
			else:
				aspect_ratio = 1000

				
		group_lengths_array.append((u_length, v_length))
		group_centres_array.append(centre)
		group_aspect_array.append(aspect_ratio)


	#find the group with the highest aspect ratio and make that the starting point
	highest_aspect = max(group_aspect_array)
	highest_aspect_index = group_aspect_array.index(highest_aspect)
	it = highest_aspect_index


	roundvaluex = 0
	roundvaluey = 0

	for i in range(len(edge_groups)):	#starting with highest aspect ratio group then continuing through them all
		group = edge_groups[it]

		u_length = group_lengths_array[it][0]
		v_length = group_lengths_array[it][1]
		centre = group_centres_array[it]


		if len(group) == 0:					#keeps the lockstep working even when there's empty groups
			if snap_direction == "U":
				snap_direction == "V"
			else:
				snap_direction == "U"


		if snap_direction == None:			#this only runs on the first group, everything else continues in lockstep from there
			if u_length < v_length :
				snap_direction = "U"
			else:
				snap_direction = "V"

		elif snap_direction == "U":			#if the snap direction has been used already we proceed in lockstep, alternating between the axes
			snap_direction = "V"
		elif snap_direction == "V":
			snap_direction = "U"


		if snap_direction == "U":
			if roundvaluex == 0 and len(group) != 0:	
				if centre[0] < UVBoundsCentre[0]:		#if the centre of the snapping group is less than the centre of the shell, round downwards, otherwise up - only for first in each axis, then lockstep
					roundvaluex = -1
				else:
					roundvaluex = 1
			elif roundvaluex == 1:					#flip them round if already used
				roundvaluex = -1
			elif roundvaluex == -1:
				roundvaluex = 1

			if roundvaluex == 1:
				snap_value = CurrentBox.right_bound		#pick which bound to snap to based on the round value
			else:
				snap_value = CurrentBox.left_bound


		if snap_direction == "V":
			if roundvaluey == 0 and len(group) != 0:	
				if centre[1] < UVBoundsCentre[1]:		#if the centre of the snapping group is less than the centre of the shell, round downwards, otherwise up - only for first in each axis, then lockstep
					roundvaluey = -1
				else:
					roundvaluey = 1
			elif roundvaluey == 1:					#flip them round if already used
				roundvaluey = -1
			elif roundvaluey == -1:
				roundvaluey = 1

			if roundvaluey == 1:
				snap_value = CurrentBox.top_bound		#now using box bounds instead of doing a comparison
			else:
				snap_value = CurrentBox.bottom_bound


		for vert in group:
			if snap_direction == "U":
				vert.co.x = snap_value	#snap vert to bounds in specified axis
			else:
				vert.co.y = snap_value


		if it < (len(edge_groups)-1):	#lets us loop back round to the start
			it += 1
		else:
			it = 0

	return bm

def CleanupWorkingDuplicate(old_bm, uvlayer, og_obj, og_uv_locs, og_uv_seams, og_selection, og_mesh_select_mode, og_uv_sync, og_uv_select_mode): #transfers UV movement back to the original mesh, deletes the duplicate

	#create a dictionary to hold the original locs of all the uv loops and their new position.
	#This is currently done purely by uv pos	    
	uv_loc_dict = {}


	#loop through all verts
	for vert in old_bm.verts:
		new_coords = vert.co.xy
		og_coords = og_uv_locs[vert]
		
		#add new coords to dict with old as key (needs to be a tuple for some reason)
		uv_loc_dict[tuple(og_coords)] = new_coords
	 
	#cleanup duplicate    
	bpy.ops.object.editmode_toggle()  
	bpy.ops.object.delete(use_global=False)
	  
	#change active object to the original    
	og_obj.select_set(True)
	bpy.context.view_layer.objects.active = og_obj


	#bmesh setup using original this time
	obj = bpy.context.view_layer.objects.active 
	me = obj.data

	#get back into edit mode
	if bpy.context.object.mode == 'OBJECT':
		bpy.ops.object.editmode_toggle() 
		
	bm = bmesh.from_edit_mesh(me) 
	uvlayer = bm.loops.layers.uv.verify()

	bpy.context.tool_settings.uv_select_mode = og_uv_select_mode

	#iterate through verts - ideally limit this to the ones we know are selected, but we can't use tagging or anything, so it would have to be some kind of "not in" array which isn't great
	for vert in bm.verts:
		for loop in vert.link_loops:
			#for every loop, try to grab the uv value, stick it into the dict and change coords to the "new value"
			try:
				oldvalue = loop[uvlayer].uv
				newvalue =  uv_loc_dict[tuple(oldvalue)]
				loop[uvlayer].uv = newvalue
				loop[uvlayer].select = True
			except:
				continue

	
	for edge in bm.edges:
		if edge.index in og_uv_seams:	#we have to use index here because it kills the original bmesh really early for some reason
			edge.seam = True


	#makes sure all the indices are accurate
	bm.faces.ensure_lookup_table()
	bm.verts.ensure_lookup_table()
	bm.edges.ensure_lookup_table()

	#restoring original selections
	bpy.ops.mesh.select_all(action='DESELECT')

	bpy.context.tool_settings.mesh_select_mode[:] = og_mesh_select_mode
	bpy.context.tool_settings.use_uv_select_sync = og_uv_sync

	if og_mesh_select_mode[0]== True:
		for index in og_selection[0]:
			vert = bm.verts[index]
			vert.select = True
	if og_mesh_select_mode[1] == True:
		for index in og_selection[1]:
			edge = bm.edges[index]
			edge.select = True
			for vert in edge.verts:
				vert.select = True
	if og_mesh_select_mode[2] == True:
		for index in og_selection[2]:
			face = bm.faces[index]
			face.select = True 

	if og_mesh_select_mode[2] == False:
		for face in bm.faces:
			face_selected = True
			for vert in face.verts:
				if vert.select == False:
					face_selected = False 
			if face_selected == True:
				face.select = True



	bmesh.update_edit_mesh(og_obj.data)

	return og_obj

def SoftSelect(bm, corners, pre_shellscale_corner_locs, UVBounds, vert_array, CurrentBox, pre_correction_corner_locs, vert_og_locs_dict, pre_snap_corners_locs, context): #moves all vertices based on how close they are to each corner and how far the corners moved
	
	#There is a lot of stuff in here. Most of it is there to ensure that soft select does something, but doesn't push any vertices outside of the box bounds

	if bpy.context.scene.use_full_soft_select == False:
		print("Skipping soft select")
		return bm

	corners_og_locs = []
	corner_trans_array= []

	#create noncorners array by removing corners from vert_array
	noncorners = vert_array
	for vert in corners:
		try:
			noncorners.remove(vert)
		except:
			continue

	for vert in corners:
		startxloc = pre_snap_corners_locs[vert.index][0] 
		startyloc = pre_snap_corners_locs[vert.index][1] 
		corners_og_locs.append([startxloc, startyloc])		#we need an array instead of a vector I think?

		corner_trans = numpy.array([((vert.co.x - startxloc)* context.scene.softselect_multiplier),((vert.co.y - startyloc)* context.scene.softselect_multiplier)])	#find the distance between the corner before and after snapping - this is the transform we use for the soft select
		corner_trans_array.append(corner_trans)


	vert_to_corner_distance = {} #vert key, array of distances value

	for vert in noncorners:
		dist_to_corners = []
		for i in range(4):
			corner = corners[i]
			og_vert_loc = vert_og_locs_dict[vert]
			true_distance = math.dist(og_vert_loc,pre_correction_corner_locs[corner])
			dist_to_corners.append(abs(true_distance)) 		
		vert_to_corner_distance[vert] = dist_to_corners


	maxdists = []	#distance values from each corner to the furthest away vert
	for i in range(4):	#finding max distance from all four corners
		array = []
		corner = corners[i]
		for dist in dist_to_corners:
			array.append(dist)
		maxdist = max(array)
		maxdists.append(maxdist * context.scene.softselect_falloff)


	vert_to_scaled_distances = {} #scaled distance is max distance for that corner divided by the vertices distance, inverted
	vert_to_axis_transforms = {}

	for vert in noncorners:
		scaled_distances = []
		for i in range(4):
			scaled_distance = Saturate(1 - numpy.linalg.norm((vert_to_corner_distance[vert][i] / maxdists[i]))) #scaled distance values for all four corners
			scaled_distances.append(scaled_distance)
		vert_to_scaled_distances[vert] = scaled_distances	#this basically means we've got a decimal value to multiply each corner transform by for each vert

		x_positive = 0
		y_positive = 0
		y_negative = 0 
		x_negative = 0

		for i in range(4):
			corner_transform = corner_trans_array[i] * vert_to_scaled_distances[vert][i]	#transform is corner transform multiplied by the scaled dist for that ocrner
			corner_transform_x = corner_transform[0]
			corner_transform_y = corner_transform[1]

			#split each transform up to axis and signs to find the total transform in each cardinal direction
			if corner_transform_x >= 0:	
				x_positive += corner_transform_x
			if corner_transform_x < 0:
				x_negative += corner_transform_x
			if corner_transform_y >= 0:		
				y_positive += corner_transform_y
			if corner_transform_y < 0:
				y_negative += corner_transform_y

		axis_transforms = [x_positive, y_positive, x_negative, y_negative]	#axis transforms for each vertex
		vert_to_axis_transforms[vert] = axis_transforms

	axis_scalars = [] #multiplied against the axis transform for each vertex

	for i in range(4):			#once for each cardinal direction/axis

		highest_oob_dist = 0
		highest_bound_direction = 0
		highest_bound_axis = None
		oob_vert = None

		for vert in noncorners:

			#find the new position of the vertex in the current axis
			if i == 0 or i == 2:
				newpos = [vert.co.x + vert_to_axis_transforms[vert][i], vert.co.y]
			if i == 1 or i ==3:
				newpos = [vert.co.x, vert.co.y + vert_to_axis_transforms[vert][i]]

			#check to see if it's out of bounds
			oob_dist, bound_value, bound_direction, bound_axis = OutOfBounds(newpos, CurrentBox)

			#find the most out of bounds vertex
			if abs(oob_dist) > abs(highest_oob_dist):
				highest_oob_dist = oob_dist
				oob_vert = vert
				highest_bound_value = bound_value
				highest_bound_direction = bound_direction
				highest_bound_axis = bound_axis

		print("\n Highest out of bounds vert found: distance =", highest_oob_dist, "direction = ",highest_bound_direction, "highest_bound_axis",highest_bound_axis)

		if oob_vert == None:	#if nothing is out of bounds we can leave it alone
			print("No out of bounds found ,scalar left at 1")
			scalar = 1

		if oob_vert != None:	#scalar values required to keep all verts in the atlas box
			oob_vert_transform = vert_to_axis_transforms[oob_vert][i]
			padding = 0 #CurrentBox.length/50
			if highest_bound_direction == 1 :
				if highest_bound_axis == "X":
					scalar  = ((highest_bound_value - padding) - oob_vert.co.x)/oob_vert_transform
				if highest_bound_axis == "Y":
					scalar  =  ((highest_bound_value - padding) - oob_vert.co.y)/oob_vert_transform
			if highest_bound_direction == -1:
				if highest_bound_axis == "X":
					scalar =((highest_bound_value + padding) - oob_vert.co.x)/ oob_vert_transform
				if highest_bound_axis == "Y":
					scalar = ((highest_bound_value + padding) - oob_vert.co.y) /oob_vert_transform
		axis_scalars.append(scalar)

	print("Axis axis_scalars", axis_scalars)

	for vert in noncorners: #apply scaled transform
		final_x_trans = 0 
		final_y_trans = 0
		for i in range(4):
			if i == 0 or i == 2:
				final_x_trans += vert_to_axis_transforms[vert][i] * axis_scalars[i]
			if i == 1 or i ==3:
				final_y_trans += vert_to_axis_transforms[vert][i] * axis_scalars[i]
		vert.co.x += final_x_trans
		vert.co.y += final_y_trans

	
	return bm
	
#~~~~~~~~~~The actual operator~~~~~~~~~~~

class SnapUVToAtlas (bpy.types.Operator):
	bl_idname = "uv.snaptoatlas"
	bl_label = "Snap UV to Atlas"
	bl_options = {"REGISTER", "UNDO"} 

	def execute(self, context):

		print("------------------------START------------------------")

		start = timer()

		#validate selection to single mesh
		meshes = []
		for obj in bpy.context.selected_objects:
			if obj.type == "MESH":
				meshes.append(obj)

		if len(meshes) > 1:
			self.report({"ERROR"}, "More than one mesh selected!")
			return {'CANCELLED'}

		#read the scene atlas and get it's uv coords as two lists of floats 
		atlas_u, atlas_v = ReadAtlas(context)

		#splits off the originally selected uvs on the actual mesh, sets up a bmesh copy of the original for all further operations
		bm, uvlayer , og_obj, og_uv_sync = SetupBmesh(context) 

		#Saves some states that will be changed later so we can reinstate them at the end
		og_mesh_select_mode, og_uv_seams, og_selection, og_uv_select_mode = SaveInitialState(bm, uvlayer)


		#Limits selected UVs to ones that are currently visible 
		bm, SelectedUVs, selected_uv_loops = FindSelectedUVs(bm, uvlayer, context)


		#Creates a duplicate of the current mesh and makes it into a topological representation of the uv selection
		bm,  og_uv_locs,   edge_verts_array, edge_uv_array, vert_array, og_edge_verts_locs, is_shell_only_corners = CreateWorkingDuplicate(bm, SelectedUVs, uvlayer, og_uv_seams)

		#Find the atlas box that the selected UVs lie within
		CurrentBox, UVBoundsCentre, UVBounds = FindCurrentBox(bm, vert_array, edge_uv_array, atlas_u, atlas_v) 


		#Finds corners by looking for closest value to box corners. Snaps them to the location of those corners. 
		vert_array, corners, pre_shellscale_corner_locs, start_corner_ind, UVBounds, pre_correction_corner_locs, vert_og_locs_dict, pre_snap_corners_locs = FindAndSnapCorners(bm, uvlayer, vert_array, CurrentBox, edge_verts_array, UVBounds, context)


		if is_shell_only_corners == False: #we only need to do this if there are more than four verts 

			starting_vert = corners[start_corner_ind] #pick the closest corner as the starting vert for consistency

			#apply soft select function to all verts excluding corners
			bm = SoftSelect(bm, corners, pre_shellscale_corner_locs, UVBounds, vert_array, CurrentBox, pre_correction_corner_locs, vert_og_locs_dict, pre_snap_corners_locs, context)

			#march from vert to vert using seamed edges
			bm = MarchEdges(bm, corners, edge_verts_array, atlas_u, atlas_v, UVBoundsCentre, starting_vert, og_edge_verts_locs, CurrentBox)


		og_obj = CleanupWorkingDuplicate(bm, uvlayer, og_obj, og_uv_locs, og_uv_seams, og_selection, og_mesh_select_mode, og_uv_sync, og_uv_select_mode)

		end = timer()
		print("Whole program took", end - start, "seconds")

		return {'FINISHED'}



#------Blender Requirements------

CLASSES = [                 #list of all classes
	UV_Snapper_UV_PT_Panel,
	UV_Snapper_Settings_UV_PT_Panel,
	SnapUVToAtlas
]

def register():
	for (prop_name, prop_value) in PROPS:
		setattr(bpy.types.Scene, prop_name, prop_value)
	
	for klass in CLASSES:
		bpy.utils.register_class(klass)

def unregister():
	for (prop_name, _) in PROPS:
		delattr(bpy.types.Scene, prop_name)

	for klass in CLASSES:
		bpy.utils.unregister_class(klass)
		

if __name__ == '__main__':
	register()



