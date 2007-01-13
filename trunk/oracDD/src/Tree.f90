!!$/---------------------------------------------------------------------\
!!$                                                                      |
!!$  Copyright (C) 2006-2007 Massimo Marchi <Massimo.Marchi@cea.fr>      |
!!$                                                                      |
!!$      This program is free software;  you  can  redistribute  it      |
!!$      and/or modify it under the terms of the GNU General Public      |
!!$      License version 2 as published  by  the  Free  Software         |
!!$      Foundation;                                                     |
!!$                                                                      |
!!$      This program is distributed in the hope that  it  will  be      |
!!$      useful, but WITHOUT ANY WARRANTY; without even the implied      |
!!$      warranty of MERCHANTABILITY or FITNESS  FOR  A  PARTICULAR      |
!!$      PURPOSE.   See  the  GNU  General  Public License for more      |
!!$      details.                                                        |
!!$                                                                      |
!!$      You should have received a copy of the GNU General  Public      |
!!$      License along with this program; if not, write to the Free      |
!!$      Software Foundation, Inc., 59  Temple  Place,  Suite  330,      |
!!$      Boston, MA  02111-1307  USA                                     |
!!$                                                                      |
!!$\---------------------------------------------------------------------/
! (c) Copyright Michael Metcalf and John Reid, 1992. This file may be
! freely used and copied for educational purposes provided this notice
! remains attached. Extracted from "Fortran 90 Explained" Oxford
! University Press (Oxford and New York), ISBN 0-19-853772-7.
!
module Tree
!
! Strong typing imposed
  USE Constants, ONLY: max_char_tree,max_char_long
  USE Strings, ONLY: MY_Fxm
  USE Errors, ONLY: abort_now
  implicit none
!
! Only subroutine interfaces, the length of the character
! component, and the I/O unit number are public
   private
   public Tree__start, Tree__add_node, Tree__remove_node, Tree__retrieve, Tree__finish&
        &, Tree__get_tree, Tree__print_tree, trees, branch, Tree__check_tree
!
! Module constants
   character(*), parameter:: eot = 'End-of-Tree.....'
!
! Define the basic tree type
   type trees
       character(max_char_tree) :: name    ! name of node
       CHARACTER(len=max_char_long), POINTER  :: y(:)    ! stored real data
       type(trees), pointer :: parent  ! parent node
       type(trees), pointer :: sibling ! next sibling node
       type(trees), pointer :: child   ! first child node
    end type trees
  TYPE branch
     CHARACTER(max_char_tree)   :: name,parent
     CHARACTER(max_char_long), DIMENSION(:), POINTER :: data
     CHARACTER(max_char_tree), DIMENSION(:), POINTER :: children=>NULL()
  END TYPE branch
!
! Module variables
   type(trees), pointer    :: current       ! current node
   type(trees), pointer    :: forest_root   ! the root of the forest
   integer                :: max_data      ! max size of data array
   character(max_char_tree), allocatable, target :: names(:)
                                           ! for returning list of names
! The module procedures
 
contains 
   subroutine Tree__start(My_Root)
     TYPE(trees), POINTER :: My_root
! Initialize the tree.
      allocate (forest_root)
      current => forest_root
      forest_root%name = 'forest_root'
      nullify(forest_root%parent, forest_root%sibling, forest_root%child)
      allocate(forest_root%y(0))
      max_data = 0
      IF(allocated(names)) THEN
         DEALLOCATE(names)
      END IF
      allocate (names(0))
      my_root=>forest_root
   end subroutine Tree__start
   subroutine Tree__find(name)
      character(*), intent(in) :: name
! Make the module variable current point to the node with given name,
! or be null if the name is not there.
      type(trees), pointer    :: root
! For efficiency, we search the tree rooted at current, and if this
! fails try its parent and so on until the forest root is reached.
      if (associated(current)) then
         root => current
         nullify (current)
      else
         root => forest_root
      end if
      do
         call look(root, name)
         if (associated(current)) return
         root => root%parent
         if (.not.associated(root)) exit
      end do
   contains
      recursive subroutine look(root, name)
         character(*), intent(in)        :: name
         type(trees), intent(in), target  :: root
! Look for name in the tree rooted at root. If found, make the
! module variable current point to the node
         type(trees), pointer    :: child
!
         if (MY_Fxm(TRIM(name),TRIM(root%name))) then
            current => root
         else
            child => root%child
            do
               if (.not.associated(child)) exit
               call look(child, name)
               if (associated(current)) return
               child => child%sibling
            end do
         end if
      end subroutine look
   end subroutine Tree__find
   subroutine Tree__add_node(name, name_of_parent, data)
      character(*), intent(in)   :: name, name_of_parent
! For a root, name = ''
    CHARACTER(len=*), intent(in), optional :: data(:)
! Allocate a new tree node of type node, store the given name and
! data there, set pointers to the parent and to its next sibling
! (if any). If the parent is not found, the new node is treated as
! a root. It is assumed that the node is not already present in the
! forest.
      type(trees), pointer :: new_node
!
      allocate (new_node)
      new_node%name = name
      if (present(data)) then
         allocate(new_node%y(size(data)))
         new_node%y = data
         max_data = max(max_data, size(data))
      else
         allocate(new_node%y(0))
      end if
!
! If name of parent is not null, search for it. If not found, print message.
      if (name_of_parent == '') then
         current => forest_root
      else
         call Tree__find (name_of_parent)
         if (.not.associated(current)) then
            print *, 'no parent ', name_of_parent, ' found for ', name
            current => forest_root
         end if
      end if
      new_node%parent => current
      new_node%sibling => current%child
      current%child => new_node
      nullify(new_node%child)
   end subroutine Tree__add_node
 
   subroutine Tree__remove_node(name)
      character(*), intent(in) :: name
! Remove node and the subtree rooted on it (if any),
! deallocating associated pointer targets.
      type(trees), pointer :: parent, child, sibling
!
      call Tree__find (name)
      if (associated(current)) then
         parent =>  current%parent
         child => parent%child
         if (.not.associated(child, current)) then
! Make it the first child, looping through the siblings to find it
! and resetting the links
            parent%child => current
            sibling => child
            do
              if (associated (child%sibling, current)) exit
              sibling => sibling%sibling
            end do
            sibling%sibling => current%sibling
            current%sibling => child
         end if
         call Tree__remove(current)
      end if
   end subroutine Tree__remove_node
   recursive subroutine Tree__remove (old_node)
! Remove a first child node and the subtree rooted on it (if any),
! deallocating associated pointer targets.
      type(trees), pointer :: old_node
      type(trees), pointer :: child, sibling
!
      child => old_node%child
      do
         if (.not.associated(child)) exit
         sibling => child%sibling
         call Tree__remove(child)
         child => sibling
      end do
! remove leaf node
      if (associated(old_node%parent)) old_node%parent%child => old_node%sibling
      deallocate (old_node%y)
      deallocate (old_node)
    end subroutine Tree__remove
 
    subroutine Tree__retrieve(name, real_name, data, parent, children)
      character(*), intent(in)         :: name
      CHARACTER(len=max_char_long), pointer :: data(:)
      character(max_char_tree), intent(out) :: parent,Real_Name
      character(max_char_tree), pointer     :: children(:)
! Returns a pointer to the data at the node, the name of the
! parent, and a pointer to the names of the children.
      integer count, i
      type(trees), pointer :: child
!
      call Tree__find (name)
      if (associated(current)) then
         Real_name=current%name
         data => current%y
         parent = current%parent%name
! count the number of children
         count = 0
         child => current%child
         do
           if (.not.associated(child)) exit
           count = count + 1
           child => child%sibling
         end do
         deallocate (names)
         allocate (names(count))
! and store their names
         children => names
         child => current%child
         do i = 1, count
            children(i) = child%name
            child => child%sibling
         end do
      else
         nullify(data)
         parent = ''
         nullify(children)
      end if
   end subroutine Tree__retrieve
  
   subroutine Tree__finish
! Deallocate all allocated targets.
      call Tree__remove (forest_root)
      deallocate(names)
    end subroutine Tree__finish
   SUBROUTINE Tree__Get_tree(My_Root)
     TYPE(trees), POINTER :: My_root
     current=>My_Root
   END SUBROUTINE Tree__Get_tree
   recursive subroutine Tree__print_tree(name)
! To print the data contained in a subtree
      character(*) :: name
      integer                             i
      CHARACTER(len=max_char_long), pointer                    :: data(:)
      character(max_char_tree)                 parent, self, real_name
      character(max_char_tree), pointer     :: children(:)
      character(max_char_tree), allocatable :: siblings(:)
!
      call Tree__retrieve(name, real_name, data, parent, children)
      if (.not.associated(data)) return
      self = real_name; write(*,*) TRIM(self), (TRIM(data(i)),i=1,SIZE(data))
      write(*,*) '   parent:   ', TRIM(parent)
      if (size(children) > 0 ) write(*,*) '   children: ', ('!',TRIM(children(i)),'!',i=1,SIZE(children))
      allocate(siblings(size(children)))
      siblings = children
      do i = 1, size(children)
         call Tree__print_tree(siblings(i))
      end do
   end subroutine Tree__print_tree
   SUBROUTINE Tree__Check_Tree(name,check)
!!$! To print the data contained in a subtree
     TYPE(Branch)  :: check
     CHARACTER(len=*) :: name
     CHARACTER(len=max_char_long), POINTER  :: data(:)
     CHARACTER(max_char_tree) :: parent, self, Real_name
     CHARACTER(max_char_tree), DIMENSION(:), POINTER :: children
!!$!
     IF(ASSOCIATED(check%children)) THEN
        DEALLOCATE(check%children)
        NULLIFY(check%children)
     END IF
     call Tree__retrieve(name, Real_name, data, parent, children)
     IF(ASSOCIATED(data)) THEN
        ALLOCATE(check%children(SIZE(children)))
        check%name=Real_Name
        check%parent=parent
        check%children=children
!!$
!!$--- Stop at children and parent
!!$
     END IF
   END SUBROUTINE Tree__Check_Tree
 end module Tree
