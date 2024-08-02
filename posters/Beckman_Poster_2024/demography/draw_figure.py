# import demesdraw
import demes
import demesdraw

colors = {'Ancestral':'crimson', 'Pre':'gold', 
          'Iran':'limegreen', 
          'France':'dodgerblue',
          'Germany':'darkorange', 
          'Heligoland':'purple'}

graph = demes.load('combined_demes_draw.yaml')

ax = demesdraw.tubes(graph, colours = colors, title='Combined Demographic Inference')
ax.figure.savefig('combined_new_colors.png')

###############
colors = {'Ancestral':'crimson', 'Pre':'gold', 
          'Iran':'limegreen', 
          'France':'dodgerblue'}

graph = demes.load('IRA_FRA_demes_draw.yaml')

ax = demesdraw.tubes(graph, colours = colors, title='Iran vs. France Demographic Inference')
ax.figure.savefig('IRA_FRA_new_colors.png')

###############
colors = {'Ancestral':'crimson', 'Pre':'gold',  
          'France':'dodgerblue',
          'Germany':'darkorange'}

graph = demes.load('FRA_GER_demes_draw.yaml')

ax = demesdraw.tubes(graph, colours = colors, title='France vs. Germany Demographic Inference')
ax.figure.savefig('FRA_GER_new_colors.png')

###############
colors = {'Ancestral':'crimson', 'Pre':'gold',
          'Germany':'darkorange', 
          'Heligoland':'purple'}

graph = demes.load('GER_HEL_demes_draw.yaml')

ax = demesdraw.tubes(graph, colours = colors, title='Germany vs. Heligoland Demographic Inference')
ax.figure.savefig('GER_HEL_new_colors.png')