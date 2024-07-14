# UV Snapper

![example_image](https://github.com/user-attachments/assets/6003ccdc-1f13-4805-a77c-c7118cd42338)

## Who is this for?
Anyone that uses **texture atlases or trim sheets** which need UV shells to fill boxes – **particularly anyone that uses DreamUV to support this workflow**

The kind of assets that use these workflows are typically hard surface, but the addon will fit any shape to the box (like Suzanne as shown above!) Straight, clean shells will generally provide the best results.

## Quick Start
 * **Download "uv_snapper.py" from this repo.** Optionally, download the demo scene and atlas texture
* In Blender, **install the downloaded file** from Preferences/Add-ons -> Install from File
* In the UV Editor, open the **UV Snapper panel**, select a shell or a group of UVs and hit the **Snap UV to Atlas button** 
> You can select UVs in any of the four selection modes - as long as you have a whole face selected it will rip your selection into a new shell and snap it. UV Sync mode is also supported.

> The UVs the operator chooses to be the shell corners are the UVs closest to the box corners - if the operator picks the wrong one you can make it clearer by moving the desired corner UVs closer to the box corners
## What does it do?
**Reads a user defined atlas and snaps the current UV selection to the edges of the atlas box that it lies within**. 

Soft select helps to create a more natural spread of UVs than just snapping. There are settings provided to tweak the soft select.

![image](https://github.com/user-attachments/assets/0d244a14-bd06-4f21-bcbc-817f12d53d2f)

## Why make this?
* When using DreamUV I found that some shell shapes would not unwrap well, resulting in me having to straighten them and manually position them so that the edges lined up with the box. This was designed to solve that problem. 

* Beyond that, it’s been a huge exercise in addon development and python. I’ve rewritten the whole thing several times as my skills have improved from experience working on this and other addons.

* Big shoutout to leukbaars' [Dream UV](https://github.com/leukbaars/DreamUV) for being the jumping off point for this whole thing.

* The entire addon is commented in detail - any and all feedback is greatly appreciated. You can find me on twitter/x [here](https://twitter.com/DRReadle)
