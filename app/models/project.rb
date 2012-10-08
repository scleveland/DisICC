# Sensors Demo
# Copyright (c) 2009 Montana State University
#
# FILE: project.rb
# The actual project
class Project
  include DataMapper::Resource
  
  property :id,           Serial
  property :name,         String, :required => false
  property :description,  String, :required => false
  property :synopsis,     Text
  
  # has n, :users, :through => Resource
  # has n, :sites
  #   has n, :sensor_values, :through => :sites
  #   
  # has n, :roles
  # has n, :users, :through => :roles

  
  validates_is_unique :name
  
  default_scope(:default).update(:order => [:name]) # set default order
  
  # This method allows us to do things like
  #    yogo_project_path(@project)
  # Instead of having to put @project.id
  def to_param
    self.name
  end
  
  def get_location
    @_glat ||= sites.inject(0){|sum,site| sum + site.lat } / sites.length
    @_glong ||= sites.inject(0){|sum,site| sum + site.long } / sites.length
    return @_glat,@_glong
  end
end
