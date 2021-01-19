#!/bin/bash

xacro --inorder hand_fingers_fixed_prismatic_xyz_revolute_z.urdf.xacro > hand_fingers_fixed_prismatic_xyz_revolute_z.urdf
xacro --inorder hand_prismatic_xyz_revolute_z.urdf.xacro > hand_prismatic_xyz_revolute_z.urdf
xacro --inorder hand.urdf.xacro > hand.urdf
xacro --inorder mug_panda_arm_hand_fingers_fixed.urdf.xacro > mug_panda_arm_hand_fingers_fixed.urdf
xacro --inorder panda_arm_hand.urdf.xacro > panda_arm_hand.urdf
xacro --inorder panda_arm_hand_fingers_fixed.urdf.xacro > panda_arm_hand_fingers_fixed.urdf
xacro --inorder panda_arm.urdf.xacro > panda_arm.urdf
xacro --inorder table_mug_panda_arm_hand_fingers_fixed.urdf.xacro > table_mug_panda_arm_hand_fingers_fixed.urdf
xacro --inorder table_panda_arm_hand_fingers_fixed.urdf.xacro > table_panda_arm_hand_fingers_fixed.urdf

cd sphere_col
xacro --inorder hand_fingers_fixed_prismatic_xyz_revolute_z.urdf.xacro > hand_fingers_fixed_prismatic_xyz_revolute_z.urdf
xacro --inorder hand_prismatic_xyz_revolute_z.urdf.xacro > hand_prismatic_xyz_revolute_z.urdf
xacro --inorder hand.urdf.xacro > hand.urdf
xacro --inorder panda_arm_hand.urdf.xacro > panda_arm_hand.urdf
xacro --inorder panda_arm_hand_fingers_fixed.urdf.xacro > panda_arm_hand_fingers_fixed.urdf
xacro --inorder panda_arm.urdf.xacro > panda_arm.urdf
xacro --inorder table_mug_panda_arm_hand_fingers_fixed.urdf.xacro > table_mug_panda_arm_hand_fingers_fixed.urdf
xacro --inorder table_panda_arm_hand_fingers_fixed.urdf.xacro > table_panda_arm_hand_fingers_fixed.urdf

cd ../sphere_vis
xacro --inorder hand_fingers_fixed_prismatic_xyz_revolute_z.urdf.xacro > hand_fingers_fixed_prismatic_xyz_revolute_z.urdf
xacro --inorder hand_prismatic_xyz_revolute_z.urdf.xacro > hand_prismatic_xyz_revolute_z.urdf
xacro --inorder hand.urdf.xacro > hand.urdf
xacro --inorder mug_hand_fingers_fixed_prismatic_xyz_revolute_z.urdf.xacro > mug_hand_fingers_fixed_prismatic_xyz_revolute_z.urdf
xacro --inorder mug_hand_prismatic_xyz_revolute_z.urdf.xacro > mug_hand_prismatic_xyz_revolute_z.urdf
xacro --inorder mug_panda_arm_hand_fingers_fixed.urdf.xacro > mug_panda_arm_hand_fingers_fixed.urdf
xacro --inorder panda_arm_hand.urdf.xacro > panda_arm_hand.urdf
xacro --inorder panda_arm_hand_fingers_fixed.urdf.xacro > panda_arm_hand_fingers_fixed.urdf
xacro --inorder panda_arm.urdf.xacro > panda_arm.urdf
xacro --inorder table_mug_panda_arm_hand_fingers_fixed.urdf.xacro > table_mug_panda_arm_hand_fingers_fixed.urdf
xacro --inorder table_panda_arm_hand_fingers_fixed.urdf.xacro > table_panda_arm_hand_fingers_fixed.urdf