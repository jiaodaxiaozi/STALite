# STALite
 Light-weight computational engine for Static Traffic Assignment on large-scale transportation network
 
 STAlite is an open-source AMS library for efficiently macroscopic traffic assignment 
 based on General Modeling Network Specification (GMNS) format
  
  ![nexta](doc/images/nexta.png)
 #Features
 ##Network representation based on GMNS format
 ##Easy to include demand from different multiple time periods( AM, MD, PM, NT or Hourly)
 ##Provide API for both C++ and Python interface 
 ##Efficiently multi-threading parallel computation and memory management, implemented in C++
 	Utilize up to 40 CPU cores, 200 GB of Memory for networks with more than 50K nodes
 ##Extendable Volume Delay Function (VDF) functions:
 ## Standard BPR function
 ##BPR_X function that can obtain dynamic travel time efficiently 
 ![STALite](doc/images/output.png)


